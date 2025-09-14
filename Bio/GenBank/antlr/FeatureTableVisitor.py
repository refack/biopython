# Generated from Bio/GenBank/FeatureTable.g4 by ANTLR 4.13.1
from antlr4 import *
if "." in __name__:
    from .FeatureTableParser import FeatureTableParser
else:
    from FeatureTableParser import FeatureTableParser

# This class defines a complete generic visitor for a parse tree produced by FeatureTableParser.

class FeatureTableVisitor(ParseTreeVisitor):

    # Visit a parse tree produced by FeatureTableParser#location_parse.
    def visitLocation_parse(self, ctx:FeatureTableParser.Location_parseContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by FeatureTableParser#anticodon_parse.
    def visitAnticodon_parse(self, ctx:FeatureTableParser.Anticodon_parseContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by FeatureTableParser#location.
    def visitLocation(self, ctx:FeatureTableParser.LocationContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by FeatureTableParser#location_term.
    def visitLocation_term(self, ctx:FeatureTableParser.Location_termContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by FeatureTableParser#simple_location.
    def visitSimple_location(self, ctx:FeatureTableParser.Simple_locationContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by FeatureTableParser#accession.
    def visitAccession(self, ctx:FeatureTableParser.AccessionContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by FeatureTableParser#fuzzy_lt.
    def visitFuzzy_lt(self, ctx:FeatureTableParser.Fuzzy_ltContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by FeatureTableParser#fuzzy_gt.
    def visitFuzzy_gt(self, ctx:FeatureTableParser.Fuzzy_gtContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by FeatureTableParser#anticodon_value.
    def visitAnticodon_value(self, ctx:FeatureTableParser.Anticodon_valueContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by FeatureTableParser#amino_acid.
    def visitAmino_acid(self, ctx:FeatureTableParser.Amino_acidContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by FeatureTableParser#sequence.
    def visitSequence(self, ctx:FeatureTableParser.SequenceContext):
        return self.visitChildren(ctx)



del FeatureTableParser