import numpy as np 
from cvxpy import *

# robot trajectory
class TrajectoryUncertainTiming:
    def __init__(self, K, L, total_time, overlap, acc, dec, delta_lane, x_i=5., vx_i=10., y_i=0., vy_i=0., dt=0.1, nb_lanes = 5):
        self.par = np.array([K, L, numExps, overlap, dt])
        self.par_actions = np.array([acc, dec, delta_lane])
        self.initial_state = np.array([x_i, vx_i, y_i, vy_i])
        self.p_vec = np.array([0.4, 0.2, 0.2, 0.1, 0.1])
        self.p_continuation = 0.8
        
    def gen_robot_action(self):
        K, L = self.par[:2]
        dt = self.par[-1]
        acc, dec, delta_lane = self.par_actions
        
        actions = np.zeros((2, K, L))
        
        actions[:,:,0] = np.zeros((2,K))
        accelerate = acc*np.ones((1,K))
        decelerate = - dec*np.ones((1,K))
        actions[:,:,1] = np.concatenate((accelerate, np.zeros((1,K))), axis=0)
        actions[:,:,2] = np.concatenate((decelerate, np.zeros((1,K))), axis=0)
        
        A = np.array([[1, dt], [0, 1.]])
        B = np.array([[0.], [dt]])
        C = np.transpose(np.array([(A**(K-1-t)).dot(B) for t in np.arange(K)])[:,:,0])
        state_left = np.array([[0.],[0.]])
        state_right = np.array([[delta_lane],[0.]])
        
        # cvxpy problem
        eps = 0.00001
        u = Variable(int(K))
        objective = Minimize(sum_squares(u))
        constraints = [C*u <= state_right + eps]
        constraints += [C*u >= state_right - eps]
        
        prob = Problem(objective, constraints)
        result = prob.solve()
        u = np.array(u.value)
        
        actions[:,:,3] = np.concatenate((np.zeros((1,K)), u.T), axis=0)
        actions[:,:,4] = -np.concatenate((np.zeros((1,K)), u.T), axis=0)
        
        return actions

    def update_state(self, state, action):
        K, L, numExps, overlap, dt = self.par
        A = np.array([[1., dt, 0., 0.], [0., 1., 0., 0.], [0., 0., 1., dt], [0., 0., 0., 1.]])
        B = np.array([[dt**2/2., 0.],[dt, 0.], [0., dt**2/2.], [0., dt]])
        new_state = A.dot(state) + B.dot(action)
        return new_state
                                         
    def sample_trajectory(self):
        p_vec, p_vec_left, p_vec_right = self.p_vec, self.p_vec_left, self.p_vec_right
        K, L, numExps, overlap, dt = self.par
        max_counter = self.lim_counter
        
        # initialize
        state_robot = np.zeros((4, K*numExps))
        state = self.initial_state
        state_robot[:,0] = state
        counter = 0
        counters = np.zeros(numExps)
        disturbances = self.gen_robot_action()
        js_noise = np.zeros(numExps)
        
        # propagate
        for exp_num in np.arange(numExps):
            # sample new disturbance
            if (counter == - max_counter):
                j_noise = np.random.choice(np.arange(L), p=p_vec_left)
            elif (counter == max_counter):
                j_noise = np.random.choice(np.arange(L), p=p_vec_right)
            else:
                j_noise = np.random.choice(np.arange(L), p=p_vec)
            # update counter
            if (j_noise == 3): # turn left
                counter -= 1
            elif (j_noise == 4): # turn right
                counter += 1
            js_noise[exp_num] = j_noise
            counters[exp_num] = counter
            
            actions = disturbances[:,:,j_noise]
            # play action
            for i in np.arange(K):
                action = actions[:,i]
                state = self.update_state(state, action)
                state_robot[:,exp_num*K+i] = state          
        
        overlap_state = np.zeros((4,overlap))
        state = self.initial_state
        overlap_state[:,0] = state 
        A = np.array([[1., dt, 0., 0.], [0., 1., 0., 0.], [0., 0., 1., dt], [0., 0., 0., 1.]])
        for i in np.arange(overlap-1):
            action = np.array([[0.],[0.]])
            state = A.dot(state)
            overlap_state[:,i+1] = state
        
        state_robot = np.concatenate((overlap_state, state_robot), axis=1)
        
        return state_robot, js_noise, counters