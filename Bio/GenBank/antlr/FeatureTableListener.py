# Generated from Bio/GenBank/FeatureTable.g4 by ANTLR 4.13.1
from antlr4 import *
if "." in __name__:
    from .FeatureTableParser import FeatureTableParser
else:
    from FeatureTableParser import FeatureTableParser

# This class defines a complete listener for a parse tree produced by FeatureTableParser.
class FeatureTableListener(ParseTreeListener):

    # Enter a parse tree produced by FeatureTableParser#location_parse.
    def enterLocation_parse(self, ctx:FeatureTableParser.Location_parseContext):
        pass

    # Exit a parse tree produced by FeatureTableParser#location_parse.
    def exitLocation_parse(self, ctx:FeatureTableParser.Location_parseContext):
        pass


    # Enter a parse tree produced by FeatureTableParser#anticodon_parse.
    def enterAnticodon_parse(self, ctx:FeatureTableParser.Anticodon_parseContext):
        pass

    # Exit a parse tree produced by FeatureTableParser#anticodon_parse.
    def exitAnticodon_parse(self, ctx:FeatureTableParser.Anticodon_parseContext):
        pass


    # Enter a parse tree produced by FeatureTableParser#location.
    def enterLocation(self, ctx:FeatureTableParser.LocationContext):
        pass

    # Exit a parse tree produced by FeatureTableParser#location.
    def exitLocation(self, ctx:FeatureTableParser.LocationContext):
        pass


    # Enter a parse tree produced by FeatureTableParser#location_term.
    def enterLocation_term(self, ctx:FeatureTableParser.Location_termContext):
        pass

    # Exit a parse tree produced by FeatureTableParser#location_term.
    def exitLocation_term(self, ctx:FeatureTableParser.Location_termContext):
        pass


    # Enter a parse tree produced by FeatureTableParser#simple_location.
    def enterSimple_location(self, ctx:FeatureTableParser.Simple_locationContext):
        pass

    # Exit a parse tree produced by FeatureTableParser#simple_location.
    def exitSimple_location(self, ctx:FeatureTableParser.Simple_locationContext):
        pass


    # Enter a parse tree produced by FeatureTableParser#accession.
    def enterAccession(self, ctx:FeatureTableParser.AccessionContext):
        pass

    # Exit a parse tree produced by FeatureTableParser#accession.
    def exitAccession(self, ctx:FeatureTableParser.AccessionContext):
        pass


    # Enter a parse tree produced by FeatureTableParser#fuzzy_lt.
    def enterFuzzy_lt(self, ctx:FeatureTableParser.Fuzzy_ltContext):
        pass

    # Exit a parse tree produced by FeatureTableParser#fuzzy_lt.
    def exitFuzzy_lt(self, ctx:FeatureTableParser.Fuzzy_ltContext):
        pass


    # Enter a parse tree produced by FeatureTableParser#fuzzy_gt.
    def enterFuzzy_gt(self, ctx:FeatureTableParser.Fuzzy_gtContext):
        pass

    # Exit a parse tree produced by FeatureTableParser#fuzzy_gt.
    def exitFuzzy_gt(self, ctx:FeatureTableParser.Fuzzy_gtContext):
        pass


    # Enter a parse tree produced by FeatureTableParser#anticodon_value.
    def enterAnticodon_value(self, ctx:FeatureTableParser.Anticodon_valueContext):
        pass

    # Exit a parse tree produced by FeatureTableParser#anticodon_value.
    def exitAnticodon_value(self, ctx:FeatureTableParser.Anticodon_valueContext):
        pass


    # Enter a parse tree produced by FeatureTableParser#amino_acid.
    def enterAmino_acid(self, ctx:FeatureTableParser.Amino_acidContext):
        pass

    # Exit a parse tree produced by FeatureTableParser#amino_acid.
    def exitAmino_acid(self, ctx:FeatureTableParser.Amino_acidContext):
        pass


    # Enter a parse tree produced by FeatureTableParser#sequence.
    def enterSequence(self, ctx:FeatureTableParser.SequenceContext):
        pass

    # Exit a parse tree produced by FeatureTableParser#sequence.
    def exitSequence(self, ctx:FeatureTableParser.SequenceContext):
        pass



del FeatureTableParser