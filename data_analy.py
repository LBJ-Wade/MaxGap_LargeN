import numpy as np
import os
import itertools

path = os.getcwd()

class Data_Input(object):

    def __init__(self, file_in, save_info=False):
        self.data = np.loadtxt(file_in)
        self.dim = 3.
        if self.dim != 3:
            print 'Data does not have correct dimensionality'
            exit()
        self.n_pts = self.data.shape[0]
        self.save_info = save_info

        self.quads = np.zeros(8, dtype=object)
        self.quads[0] = np.array([True, True, True])
        self.quads[1] = np.array([True, True, False])
        self.quads[2] = np.array([True, False, True])
        self.quads[3] = np.array([True, False, False])
        self.quads[4] = np.array([False, True, True])
        self.quads[5] = np.array([False, True, False])
        self.quads[6] = np.array([False, False, True])
        self.quads[7] = np.array([False, False, False])
        
        if save_info:
            self.output_bx_file = path + '/outputs/box_file.dat'
        return

    def solve_coords(self, point):
        vol, side = self.find_box(point)
        if self.save_info:
            np.savetxt(self.output_bx_file, np.column_stack((vol, side)))
        return side[-1]

    def distance_scan(self, point):
        distances = np.zeros((self.n_pts, 3), dtype=object)
        for i in range(self.n_pts):
            quad = self.quadrant_slvr(self.data[i], point)
            distances[i] = [self.data[i], np.sqrt(np.sum((self.data[i] - point)**2.)), quad]
        distances = distances[np.argsort(distances[:,1])]
        return distances

    def quadrant_slvr(self, data_pt, sample_pt):
        for i in range(8):
            if np.all(((data_pt - sample_pt) > 0) == self.quads[i]):
                val = i
                break
        return val

    def relvnce_tst(self, point, data_svd, data_ck):
        dist1 = np.abs(point - data_svd)
        dist2 = np.abs(point - data_ck)
        if np.all(dist1 < dist2):
            return False
        else:
            return True

    def def_quad_bndry(self, point, distances):
        quad_bndry = np.zeros(8, dtype=object)
        total_pts = []
        for i in range(8):
            quad_bndry[i] = []
        
        for i in range(self.n_pts):
            q_idx = distances[i][-1]
            if len(quad_bndry[q_idx]) == 0:
                quad_bndry[q_idx].append(distances[i][0])
                total_pts.append(distances[i][0])
            else:
                chk = True
                for j in range(len(quad_bndry[q_idx])):
                    chk = self.relvnce_tst(point, quad_bndry[q_idx][j], distances[i][0])
                    if chk == False:
                        break
                if chk == True:
                    quad_bndry[q_idx].append(distances[i][0])
                    total_pts.append(distances[i][0])

        print 'Relevant Bndry Points: {:.0f}'.format(len(total_pts))
        return quad_bndry, total_pts

    def box_calculate(self, sides):
        vol = np.abs((sides[0] - sides[3]) * (sides[1] - sides[4]) * (sides[2] - sides[5]))
        return vol

    def find_box(self, point):
        dist_vec = self.distance_scan(point)
        quads, total_vec = self.def_quad_bndry(point, dist_vec)
#        sides = np.zeros(6, dtype=object)
        combined_sides = []
        volumes = []
        for i in range(8):
            quads[i] = np.asarray(quads[i])
       
        cand_0 = np.vstack((quads[0], quads[1], quads[2], quads[3]))
        cand_0 = np.append(cand_0, [[1, point[1], point[2]]], axis=0)
        for x0 in cand_0:
            s0 = x0[0]
            cand_1 = np.vstack((quads[0], quads[1], quads[4], quads[5]))
            cand_1 = np.append(cand_1[cand_1[:, 0] < s0], [[point[0], 1, point[2]]], axis=0)
            s1 = cand_1[np.argmin(np.abs(x0[1] - cand_1[:, 1]))][1]
            
            cand_2 = np.vstack((quads[0], quads[2], quads[4], quads[6]))
            cand_2 = np.append(cand_2[(cand_2[:, 0] < s0)&(cand_2[:, 1] < s1)], [[point[0], point[1], 1]], axis=0)
            s2 = cand_2[np.argmin(np.abs(x0[2] - cand_2[:, 2]))][2]

            cand_3 = np.vstack((quads[4], quads[5], quads[6], quads[7]))
            cand_3 = np.append(cand_3[(cand_3[:, 0] < s0)&(cand_3[:, 1] < s1)&(cand_3[:, 2] < s2)], [[0., point[1], point[2]]], axis=0)
            s3 = cand_3[np.argmin(np.abs(x0[0] - cand_3[:, 0]))][0]
            
            cand_4 = np.vstack((quads[2], quads[3], quads[6], quads[7]))
            cand_4 = np.append(cand_4[(cand_4[:, 0] < s0)&(cand_4[:, 1] < s1)&(cand_4[:, 2] < s2)&(cand_4[:, 0] > s3)], [[point[0], 0., point[2]]], axis=0)
            s4 = cand_4[np.argmin(np.abs(x0[1] - cand_4[:, 1]))][1]
            
            cand_5 = np.vstack((quads[1], quads[3], quads[5], quads[7]))
            cand_5 = np.append(cand_5[(cand_5[:, 0] < s0)&(cand_5[:, 1] < s1)&(cand_5[:, 2] < s2)&(cand_5[:, 0] > s3)&(cand_5[:, 1] > s4)], [[point[0], point[1], 0.]], axis=0)
            s5 = cand_5[np.argmin(np.abs(x0[2] - cand_5[:, 2]))][2]
            
            combined_sides.append([s0, s1, s2, s3, s4, s5])
            volumes.append(self.box_calculate(combined_sides[-1]))
#        sides[0] = np.concatenate((quads[0][:, 0], quads[1][:, 0], quads[2][:, 0], quads[3][:, 0]))
#        sides[1] = np.concatenate((quads[0][:, 1], quads[1][:, 1], quads[4][:, 1], quads[5][:, 1]))
#        sides[2] = np.concatenate((quads[0][:, 2], quads[2][:, 2], quads[4][:, 2], quads[6][:, 2]))
#        sides[3] = np.concatenate((quads[4][:, 0], quads[5][:, 0], quads[6][:, 0], quads[7][:, 0]))
#        sides[4] = np.concatenate((quads[2][:, 1], quads[3][:, 1], quads[6][:, 1], quads[7][:, 1]))
#        sides[5] = np.concatenate((quads[1][:, 2], quads[3][:, 2], quads[5][:, 2], quads[7][:, 2]))
#        for i in range(3):
#            sides[i] = np.append(sides[i], 1.)
#            sides[i+3] = np.append(sides[i+3], 0.)
#
#        reduced_sds = np.zeros(5, dtype=object)
#        for i in range(len(sides[0])):
#            reduced_sds[0] = sides[1]

        #combined_sides = list(itertools.product(*sides))
#        volumes = []
#        valid_sides = []
#        print 'Combinations to Check: ', len(combined_sides)
#        for i,sde in enumerate(combined_sides):
#            if self.point_inside(sde, np.asarray(total_vec)):
#                continue
#            volumes.append(self.box_calculate(sde))
#            valid_sides.append(np.asarray(sde))
        volumes = np.asarray(volumes)
        valid_sides = np.asarray(combined_sides)
        argS = np.argsort(volumes)
        volumes = volumes[argS]
        valid_sides = valid_sides[argS]
        print 'Largest Box Found.\n Volume: {:.3e}, Sides: '.format(volumes[-1]), valid_sides[-1]
        return volumes, valid_sides

    def point_inside(self, sides, points):
        val = np.zeros(3)
        val[0] = np.sum((sides[3] < points[:, 0])&(points[:,0] < sides[0]))
        val[1] = np.sum((sides[4] < points[:, 1])&(points[:,1] < sides[1]))
        val[2] = np.sum((sides[5] < points[:, 2])&(points[:,2] < sides[2]))
        if np.all(val > 0):
            return True
        return False


