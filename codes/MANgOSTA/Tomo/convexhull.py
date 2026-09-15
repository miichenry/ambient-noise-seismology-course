
import numpy as np

class ConvexHull:

    def __init__(self, points):

        n_pts = points.shape[1]
        assert (n_pts > 5)
        centre = points.mean(1)

        angles = np.apply_along_axis(self.angle_to_point, 0, points, centre)
        pts_ord = points[:, angles.argsort()]
        pts = [x[0] for x in zip(pts_ord.transpose())]
        prev_pts = len(pts) + 1
        k = 0
        while prev_pts > n_pts:
            prev_pts = n_pts
            n_pts = len(pts)
            i = -2
            while i < (n_pts - 2):
                Aij = self.area_of_triangle(centre, pts[i], pts[(i + 1) % n_pts])
                Ajk = self.area_of_triangle(centre, pts[(i + 1) % n_pts],pts[(i + 2) % n_pts])
                Aik = self.area_of_triangle(centre, pts[i], pts[(i + 2) % n_pts])
                if Aij + Ajk < Aik: del pts[i + 1]
                i += 1
                n_pts = len(pts)
            k += 1

        # remove duplicate points
        newpts=[]
        newpts.append(pts[0])
        k=0
        for i in range(1,len(pts)):
            if pts[i][0]!=newpts[k][0] and pts[i][1]!=newpts[k][1]:
                newpts.append(pts[i])
                k=k+1

        self.ch=np.asarray(newpts)

    def angle_to_point(self, point, centre):
        delta = point - centre
        res = np.arctan(delta[1] / delta[0])
        if delta[0] < 0:
            res += np.pi
        return res

    def area_of_triangle(self, p1, p2, p3):
        return np.linalg.norm(np.cross((p2 - p1), (p3 - p1)))/2.
