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

    def set(self, coords, value, index):
        self.latest = coords
        while index > len(self.values) - 1:
            self.add(self.latest, 0)
        self.values[index] = value

    def compatible(self, coords):
        return Util.distance_eu(coords, self.latest) < self.comp_radius

    def to_string(self):
        result = '('+str(round(self.original[0], 2))+', '+str(round(self.original[1], 2)) + ')\t |\t '
        for value in self.values:
            result += str(round(value, 2)) + '\t '
        return result

    def to_dict(self):
        return {'original': self.original, 'values': self.values}


class PointDataList:

    def __init__(self, comp_radius):
        self.points = []
        self.comp_radius = comp_radius

    def pad(self):
        maxlen = 0
        for point in self.points:
            if maxlen < len(point.values):
                maxlen = len(point.values)
        for point in self.points:
            while len(point.values) < maxlen:
                point.add(point.latest, 0)

    def get_point(self, coords):
        for point in self.points:
            if point.compatible(coords):
                return point
        return None

    def set(self, coords, value, index):
        found = False
        for point in self.points:
            if point.compatible(coords) and not found:
                point.set(coords, value, index)
                found = True
        if not found:
            point = PointData(coords, self.comp_radius)
            point.set(coords, value, index)
            self.points.append(point)
        self.pad()

    def project_list(self, index):
        return [point.values[index] for point in self.points]

    def to_string(self):
        result = ''
        for point in self.points:
            result += point.to_string() + '\n'
        return result

    def to_list(self):
        result = []
        for point in self.points:
            result.append(point.to_dict())
        return result
