#!/usr/bin/env python

import rospy

import vtd_interface.msg as msg
import numpy as np
from time import sleep
from Trajectory import *
from scipy.io import savemat

def sample_traj(K, L, numExps, delta_lane, v_r):
    traj = Trajectory(K, L, int(numExps) delta_lane, x_i = 0.0, vx_i = v_r)
    dt = 0.1
    state_robot, js_noise, counters = traj.sample_trajectory()

    x = state_robot[0,:]
    y = state_robot[2,:]
    t = np.array([dt*i for i in np.arange(numExps*K)])

    robot_data_10Hz = {}
    robot_data_10Hz["state_robot"] = state_robot
    robot_data_10Hz["js_noise"] = js_noise
    robot_data_10Hz["counters"] = counters

    savemat("/home/aslsim/RiskIJRR18/robot_data_10Hz_p1_test", robot_data_10Hz)

    return x, y, t

def talker():
    pub_human = rospy.Publisher('/planner_path', msg.planner_result, queue_size=10)
    pub_scenario = rospy.Publisher('/scenario_setup', msg.scenario_start, queue_size=10)

    rospy.init_node('talker', anonymous=True)

    xr = -150.0
    yr = 15.0
    vr = 40.0
    scenario = msg.scenario_start(-160.0, yr, vr, xr, yr, vr)
    pub_scenario.publish(scenario)


    sleep(0.5)

    xs,ys,ts = sample_traj(10, 5, 60, 3.0, vr)

    robot_traj = msg.planner_result(ts,xs+xr,-ys+yr)

    pub_human.publish(robot_traj)
    sleep(0.2)
    pub_scenario.publish(scenario)

if __name__ == '__main__':
    try:
        talker()
    except rospy.ROSInterruptException:
        pass