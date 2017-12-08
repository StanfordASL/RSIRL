import numpy as np
from numpy.linalg import matrix_power as mp
from cvxpy import *

# robot trajectory
class Trajectory:
    def __init__(self, K, L, numExps, delta_lane, x_i=0., vx_i=40., y_i=2., vy_i=0., ay_i=0., dt=0.1):
        self.par = np.array([K, L, numExps, overlap, dt])
        self.par_actions = np.array([acc, dec, delta_lane])
        self.initial_state = np.array([x_i, vx_i, y_i, vy_i, ay_i])
        self.p_vec = np.array([0.3, 0.3, 0.3, 0.1])

    def gen_robot_action(self):
        K, L = self.par[:2]
        K, L = int(K), int(L)
        dt = self.par[-1]
        acc, dec, delta_lane = self.par_actions

        actions = np.zeros((2, K, L+1))

        actions[:,:,0] = np.zeros((2,K))
        accelerate = np.array([ 0.5,1,1.5,1.,0.5,0.5,0,0,0,0]).reshape((1,K))
        decelerate = np.array([-3,-6,-6,-4,-3,-2,0,0,0,0]).reshape((1,K))
        actions[:,:,1] = np.concatenate((accelerate, np.zeros((1,K))), axis=0)
        actions[:,:,2] = np.concatenate((decelerate, np.zeros((1,K))), axis=0)

        A = np.array([[1, dt, (dt**2)/2], [0, 1., dt], [0., 0., 1]])
        B = np.array([[(dt**3)/6.], [(dt**2)/2], [dt]])
        C = np.copy(B)
        for t in range(K-1):
            tt = t+1
            C = np.concatenate((mp(A, tt).dot(B), C), axis=1)
        state_left = np.array([[0.],[0.], [0.]])
        state_right = np.array([[delta_lane],[0.], [0.]])

        # cvxpy problem
        eps = 0.00001
        u = Variable(K)
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
        K, L, numExps, overlap = int(K), int(L), int(numExps), int(overlap)
        A = np.array([[1., dt, 0., 0., 0.], [0., 1., 0., 0., 0.], [0., 0., 1., dt, (dt**2)/2], [0., 0., 0., 1., dt], [0., 0., 0., 0., 1.]])
        B = np.array([[(dt**2)/2., 0.],[dt, 0.], [0., (dt**3)/6.], [0., (dt**2)/2], [0, dt]])
        new_state = A.dot(state) + B.dot(action)
        return new_state

    def sample_trajectory(self):
        p_vec = self.p_vec
        K, L, numExps, overlap, dt = self.par
        K, L, numExps, overlap = int(K), int(L), int(numExps), int(overlap)

        # robot takes random maneuvers every K time steps
        state = self.initial_state
        state_robot = np.zeros((5, K*numExps))
        state_robot[:,0] = state
        disturbances = self.gen_robot_action()
        js_noise = np.zeros(numExps)
        counters = np.zeros(numExps)

        #v_limit
        v_limit = 25.

        counter = 1 # robot starts on right lane

        # propagate trajectory
        for exp_num in np.arange(numExps):
            counters[exp_num] = counter
            j_noise = np.random.choice(np.arange(L), p=p_vec)

            if (state[1] <= v_limit): # if leader too slow, increase its speed
                #p_dec = max(0,p_vec[2]*(1.0 - (v_limit - state[1])/10.))
                #p_acc = p_vec[1] + (p_vec[2]-p_dec)
                #p_vec_red = np.array([p_vec[0], p_acc, p_dec, p_vec[3], p_vec[4]])
                #j_noise = np.random.choice(np.arange(L), p=p_vec_red)
                j_noise = 2
            elif (j_noise == 3)
                if counter == 1
                    j_noise = 4
                counter = - counter

            js_noise[exp_num] = j_noise
            actions = disturbances[:,:,j_noise]

            # play action
            for i in np.arange(K):
                action = actions[:,i]
                state = self.update_state(state, action)
                state_robot[:,exp_num*K+i] = state
        return state_robot, js_noise, counters





