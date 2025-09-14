# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Python parser for GenBank feature tables using ANTLR."""

from antlr4 import InputStream, CommonTokenStream
from .FeatureTableLexer import FeatureTableLexer
from .FeatureTableParser import FeatureTableParser
from .FeatureTableVisitor import FeatureTableVisitor

class LocationVisitor(FeatureTableVisitor):
    def visitLocation_parse(self, ctx):
        return self.visit(ctx.location())

    def visitLocation_term(self, ctx):
        if ctx.simple_location():
            return self.visit(ctx.simple_location())
        elif ctx.getChild(0).getText() == 'complement':
            return ('complement', self.visit(ctx.location_term(0)))
        elif ctx.getChild(0).getText() == 'join':
            locations = [self.visit(loc) for loc in ctx.location_term()]
            return ('join', locations)
        elif ctx.getChild(0).getText() == 'order':
            locations = [self.visit(loc) for loc in ctx.location_term()]
            return ('order', locations)

    def visitSimple_location(self, ctx):
        if ctx.DOTDOT():
            start = int(ctx.NUMBER(0).getText())
            end = int(ctx.NUMBER(1).getText())
            fuzzy_start = ctx.fuzzy_lt() is not None
            fuzzy_end = ctx.fuzzy_gt() is not None
            return ('range', start, end, fuzzy_start, fuzzy_end)
        elif ctx.CARET():
            start = int(ctx.NUMBER(0).getText())
            end = int(ctx.NUMBER(1).getText())
            return ('between', start, end)
        else:
            return ('point', int(ctx.NUMBER(0).getText()))

    def visitLocation(self, ctx):
        location_term = self.visit(ctx.location_term())
        if ctx.accession():
            accession = ctx.accession().getText()
            return ('remote', accession, location_term)
        else:
            return location_term

def parse_location(location_string):
    """Parse a GenBank location string using ANTLR."""
    input_stream = InputStream(location_string)
    lexer = FeatureTableLexer(input_stream)
    stream = CommonTokenStream(lexer)
    parser = FeatureTableParser(stream)
    tree = parser.location_parse()

    visitor = LocationVisitor()
    return visitor.visit(tree)
