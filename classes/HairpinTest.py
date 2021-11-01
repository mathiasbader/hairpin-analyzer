# ---------------------------------------------------
# --                                               --
# --    Hairpin Analyzer v0.3                      --
# --                                               --
# --    --> class HairpinTestSuite                 --
# --                                               --
# --                                               --
# --  2011-Aug-01 to 2011-Aug-23                   --
# --  by Mathias Bader (mail@mathiasbader.de)      --
# --  Universitaet des Saarlandes                  --
# ---------------------------------------------------
# 

# 
# This is a Test Suite for the Hairpin Analyzer. If
# you change any code of the software, make sure to
# run the Test Suite afterwards and check whether
# all unit test still run through, to ensure
# functionality.
#


import unittest
       

# import classes, configuration and specification
import imp
HairpinClasses = imp.load_source('HairpinData', 'HairpinClasses.py')
Configuration = imp.load_source('Configuration', '../configuration/configuration.py')
Specification = imp.load_source('Specification', 'HairpinSpec.py')

# output some general information about the software
line_software_inf_1 = "   Hairpin Analyzer v" + Specification.software_version + "   "
line_software_inf_2 = "   Test Suite"
line_software_outline = ""
for i in range(len(line_software_inf_1)):
    line_software_outline += "-"
print
print line_software_outline
print line_software_inf_1
print line_software_inf_2
print line_software_outline
print


class HairpinTestCase(unittest.TestCase):
    """Superclass of Test cases for Hairpin Analyzer"""
    def setUp(self):
        pass
    def tearDown(self):
        pass

class TestMapHairpinPositions(HairpinTestCase):
    """Test cases for map_hairpin_positions()"""
    # a number of valid mappings of hairpin data
    valid_mappings = [("00001x", "530"),
                      ("11101x", "641"),
                      ("xxx01x", "985"),
                      ("100001", "400"),
                      ("0000",   "00"),
                      ("101010", "131"),
                      ("010101", "313")
                     ]
    
    def test_even_input_count(self):
        """even position count should abort"""
        self.assertRaises(HairpinClasses.InvalidInputException, HairpinClasses.HairpinFunctionality.map_hairpin_positions, ("0011xx0x1"))
    def test_invalid_input_char(self):
        """invalid characters should abort"""
        self.assertRaises(HairpinClasses.InvalidInputException, HairpinClasses.HairpinFunctionality.map_hairpin_positions, ("0041xx0x"))
    def test_empty_input(self):
        """empty string should abort"""
        self.assertRaises(HairpinClasses.InvalidInputException, HairpinClasses.HairpinFunctionality.map_hairpin_positions, (""))
    def test_correct_mapping_examples(self):
        """try example mappings"""
        for x in self.valid_mappings:
            self.assertEqual(HairpinClasses.HairpinFunctionality.map_hairpin_positions(x[0]), x[1])

class TestDefineStructure(HairpinTestCase):
    """Test cases for define_structure()"""
    # a number of known structures
    known_structures = [("00000", "no_meth"),
                        ("10000", "only_one"),
                        ("00100", "only_one"),
                        ("00001", "only_one"),
                        ("30000", "only_one"),
                        ("00300", "only_one"),
                        ("00003", "only_one"),
                        ("40000", "only_one"),
                        ("00400", "only_one"),
                        ("00004", "only_one"),
                        ("10300", "only_one"),
                        ("30100", "only_one"),
                        ("00103", "only_one"),
                        ("00301", "only_one"),
                        ("13046", "mosaic"),
                        ("10001", "mosaic"),
                        ("10111", "mosaic"),
                        ("30003", "mosaic"),
                        ("40004", "mosaic"),
                        ("40444", "mosaic"),
                        ("10340", "mosaic"),
                        ("30140", "mosaic"),
                        ("13131", "mosaic"),
                        ("40103", "mosaic"),
                        ("11000", "continuous"),
                        ("00110", "continuous"),
                        ("00011", "continuous"),
                        ("33000", "continuous"),
                        ("00330", "continuous"),
                        ("00033", "continuous"),
                        ("44000", "continuous"),
                        ("00440", "continuous"),
                        ("00044", "continuous"),
                        ("11111", "continuous"),
                        ("33333", "continuous"),
                        ("44444", "continuous"),
                        ("14000", "continuous"),
                        ("34000", "continuous"),
                        ("00143", "continuous"),
                        ("33100", "continuous"),
                        ("41100", "continuous"),
                        ("01433", "continuous")
                       ]
    def test_define_structures(self):
        """test for definition of structures"""
        for x in self.known_structures:
            self.assertEqual(HairpinClasses.HairpinFunctionality.define_structure(x[0]), x[1])

# build the unit test suite
suite1 = unittest.TestLoader().loadTestsFromTestCase(TestMapHairpinPositions)
suite2 = unittest.TestLoader().loadTestsFromTestCase(TestDefineStructure)
alltests = unittest.TestSuite((suite1, suite2))

# run the unit test suite
unittest.TextTestRunner(verbosity=2).run(alltests)

