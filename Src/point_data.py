from util import Util


class PointData:

    def __init__(self, coords, comp_radius):
        self.original = coords
        self.latest = coords
        self.comp_radius = comp_radius
        self.values = []

    def add(self, coords, value):
        self.latest = coords
        self.values.append(value)

    def compatible(self, coords):
        return Util.distance_eu(coords, self.latest) < self.comp_radius


class PointDataList:

    def __init__(self, comp_radius):
        self.points = []
        self.comp_radius = comp_radius

    def add(self, coords, value):
        found = False
        for point in self.points:
            if point.compatible(coords) and not found:
                point.add(coords, value)
                found = True
        if not found:
            point = PointData(coords, self.comp_radius)
            point.add(coords, value)
            self.points.append(point)
