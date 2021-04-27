# coding=utf-8

class ElementData(object):
    """Record element indexes, names and masses.

    Including common elements: H-Kr, Ag, I, Au
    """

    def __init__(self):
        self.elements = {}
        self.element_indexes = {}

        # Add element infomation from periodic table
        self.add_element(1,  'H',    1.008)  # 1.00794
        self.add_element(2,  'HE',   4.003)  # 4.0026
        self.add_element(3,  'LI',   6.941)
        self.add_element(4,  'BE',   9.012)  # 9.01218
        self.add_element(5,  'B',   10.811)
        self.add_element(6,  'C',   12.011)  # 12.0107
        self.add_element(7,  'N',   14.007)  # 14.0067
        self.add_element(8,  'O',   15.999)  # 15.9994
        self.add_element(9,  'F',   18.998)  # 18.9984
        self.add_element(10, 'Ne',  20.180)  # 20.1797
        self.add_element(11, 'Na',  22.990)  # 22.9898
        self.add_element(12, 'Mg',  24.305)
        self.add_element(13, 'Al',  26.982)  # 26.9815
        self.add_element(14, 'Si',  28.086)  # 28.0855
        self.add_element(15, 'P',   30.974)  # 30.9738
        self.add_element(16, 'S',   32.065)
        self.add_element(17, 'Cl',  79.904)
        self.add_element(18, 'Ar',  39.948)
        self.add_element(19, 'K',   39.098)  # 39.0983
        self.add_element(20, 'Ca',  40.078)
        self.add_element(21, 'Sc',  44.956)  # 44.9559
        self.add_element(22, 'Ti',  47.867)
        self.add_element(23, 'V',   50.942)  # 50.9415
        self.add_element(24, 'Cr',  51.996)  # 51.9961
        self.add_element(25, 'Mn',  54.938)
        self.add_element(26, 'Fe',  55.845)
        self.add_element(27, 'Co',  58.933)  # 58.9332
        self.add_element(28, 'Ni',  58.693)  # 58.6934
        self.add_element(29, 'Cu',  63.546)
        self.add_element(30, 'Zn',  65.409)
        self.add_element(31, 'Ga',  69.723)
        self.add_element(32, 'Ge',  72.640)
        self.add_element(33, 'As',  74.922)  # 74.9216
        self.add_element(34, 'Se',  78.960)
        self.add_element(35, 'Br',  79.904)
        self.add_element(36, 'Kr',  83.798)
        self.add_element(47, 'Ag', 107.868)  # 107.8682
        self.add_element(53, 'I',  126.904)  # 126.90447
        self.add_element(79, 'Au', 196.966)  # 196.96655

    def add_element(self, index, name, mass):
        """Add element information to records."""
        self.elements[name] = {"index": index, "mass": mass}
        self.element_indexes[index] = {"name": name, "mass": mass}

    def index2element(self, i):
        """Find element name by its index."""
        if i in self.element_indexes:
            return self.element_indexes[i]["name"]
        else:
            return ''

    def index2mass(self, i):
        """Find element mass by its index."""
        if i in self.element_indexes:
            return self.element_indexes[i]["mass"]
        else:
            return 0.

    def element2index(self, element):
        """Find element name by its name."""
        if element in self.elements:
            return self.elements[element]["index"]
        else:
            return 0

    def element2mass(self, element):
        """Find element mass by its name."""
        if element in self.elements:
            return self.elements[element]["mass"]
        else:
            return 0.
