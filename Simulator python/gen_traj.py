from Trajectory import * 

def sample_traj(K, L, numExps, overlap, acc, dec, delta_lane):
	traj = Trajectory(K, L, numExps, overlap, acc, dec, delta_lane)
	dt = 0.1
	state_robot, js_noise, counters = traj.sample_trajectory()

	xs = state_robot[0,:]
	ys = state_robot[2,:]
	ts = np.array([dt*i for i in np.arange(overlap+numExps*K)])

	return xs, ys, ts

