# Generated from Bio/GenBank/FeatureTable.g4 by ANTLR 4.13.1
# encoding: utf-8
from antlr4 import *
from io import StringIO
import sys
if sys.version_info[1] > 5:
	from typing import TextIO
else:
	from typing.io import TextIO

def serializedATN():
    return [
        4,1,18,98,2,0,7,0,2,1,7,1,2,2,7,2,2,3,7,3,2,4,7,4,2,5,7,5,2,6,7,
        6,2,7,7,7,2,8,7,8,2,9,7,9,2,10,7,10,1,0,1,0,1,0,1,1,1,1,1,1,1,2,
        1,2,1,2,3,2,32,8,2,1,2,1,2,1,3,1,3,1,3,1,3,1,3,1,3,1,3,1,3,1,3,1,
        3,5,3,46,8,3,10,3,12,3,49,9,3,1,3,1,3,1,3,3,3,54,8,3,1,4,3,4,57,
        8,4,1,4,1,4,1,4,3,4,62,8,4,1,4,1,4,1,4,1,4,1,4,3,4,69,8,4,1,5,1,
        5,1,5,3,5,74,8,5,1,6,1,6,1,7,1,7,1,8,1,8,1,8,1,8,1,8,1,8,1,8,1,8,
        1,8,1,8,1,8,1,8,1,8,1,8,1,9,1,9,1,10,1,10,1,10,0,0,11,0,2,4,6,8,
        10,12,14,16,18,20,0,1,1,0,3,4,95,0,22,1,0,0,0,2,25,1,0,0,0,4,31,
        1,0,0,0,6,53,1,0,0,0,8,68,1,0,0,0,10,70,1,0,0,0,12,75,1,0,0,0,14,
        77,1,0,0,0,16,79,1,0,0,0,18,93,1,0,0,0,20,95,1,0,0,0,22,23,3,4,2,
        0,23,24,5,0,0,1,24,1,1,0,0,0,25,26,3,16,8,0,26,27,5,0,0,1,27,3,1,
        0,0,0,28,29,3,10,5,0,29,30,5,11,0,0,30,32,1,0,0,0,31,28,1,0,0,0,
        31,32,1,0,0,0,32,33,1,0,0,0,33,34,3,6,3,0,34,5,1,0,0,0,35,36,5,2,
        0,0,36,37,5,8,0,0,37,38,3,6,3,0,38,39,5,9,0,0,39,54,1,0,0,0,40,41,
        7,0,0,0,41,42,5,8,0,0,42,47,3,6,3,0,43,44,5,10,0,0,44,46,3,6,3,0,
        45,43,1,0,0,0,46,49,1,0,0,0,47,45,1,0,0,0,47,48,1,0,0,0,48,50,1,
        0,0,0,49,47,1,0,0,0,50,51,5,9,0,0,51,54,1,0,0,0,52,54,3,8,4,0,53,
        35,1,0,0,0,53,40,1,0,0,0,53,52,1,0,0,0,54,7,1,0,0,0,55,57,3,12,6,
        0,56,55,1,0,0,0,56,57,1,0,0,0,57,58,1,0,0,0,58,59,5,17,0,0,59,61,
        5,12,0,0,60,62,3,14,7,0,61,60,1,0,0,0,61,62,1,0,0,0,62,63,1,0,0,
        0,63,69,5,17,0,0,64,65,5,17,0,0,65,66,5,13,0,0,66,69,5,17,0,0,67,
        69,5,17,0,0,68,56,1,0,0,0,68,64,1,0,0,0,68,67,1,0,0,0,69,9,1,0,0,
        0,70,73,5,16,0,0,71,72,5,1,0,0,72,74,5,17,0,0,73,71,1,0,0,0,73,74,
        1,0,0,0,74,11,1,0,0,0,75,76,5,14,0,0,76,13,1,0,0,0,77,78,5,15,0,
        0,78,15,1,0,0,0,79,80,5,8,0,0,80,81,5,5,0,0,81,82,5,11,0,0,82,83,
        3,4,2,0,83,84,5,10,0,0,84,85,5,6,0,0,85,86,5,11,0,0,86,87,3,18,9,
        0,87,88,5,10,0,0,88,89,5,7,0,0,89,90,5,11,0,0,90,91,3,20,10,0,91,
        92,5,9,0,0,92,17,1,0,0,0,93,94,5,16,0,0,94,19,1,0,0,0,95,96,5,16,
        0,0,96,21,1,0,0,0,7,31,47,53,56,61,68,73
    ]

class FeatureTableParser ( Parser ):

    grammarFileName = "FeatureTable.g4"

    atn = ATNDeserializer().deserialize(serializedATN())

    decisionsToDFA = [ DFA(ds, i) for i, ds in enumerate(atn.decisionToState) ]

    sharedContextCache = PredictionContextCache()

    literalNames = [ "<INVALID>", "'.'", "'complement'", "'join'", "'order'", 
                     "'pos'", "'aa'", "'seq'", "'('", "')'", "','", "':'", 
                     "'..'", "'^'", "'<'", "'>'" ]

    symbolicNames = [ "<INVALID>", "<INVALID>", "COMPLEMENT", "JOIN", "ORDER", 
                      "POS", "AA", "SEQ", "LPAREN", "RPAREN", "COMMA", "COLON", 
                      "DOTDOT", "CARET", "LT", "GT", "IDENTIFIER", "NUMBER", 
                      "WS" ]

    RULE_location_parse = 0
    RULE_anticodon_parse = 1
    RULE_location = 2
    RULE_location_term = 3
    RULE_simple_location = 4
    RULE_accession = 5
    RULE_fuzzy_lt = 6
    RULE_fuzzy_gt = 7
    RULE_anticodon_value = 8
    RULE_amino_acid = 9
    RULE_sequence = 10

    ruleNames =  [ "location_parse", "anticodon_parse", "location", "location_term", 
                   "simple_location", "accession", "fuzzy_lt", "fuzzy_gt", 
                   "anticodon_value", "amino_acid", "sequence" ]

    EOF = Token.EOF
    T__0=1
    COMPLEMENT=2
    JOIN=3
    ORDER=4
    POS=5
    AA=6
    SEQ=7
    LPAREN=8
    RPAREN=9
    COMMA=10
    COLON=11
    DOTDOT=12
    CARET=13
    LT=14
    GT=15
    IDENTIFIER=16
    NUMBER=17
    WS=18

    def __init__(self, input:TokenStream, output:TextIO = sys.stdout):
        super().__init__(input, output)
        self.checkVersion("4.13.1")
        self._interp = ParserATNSimulator(self, self.atn, self.decisionsToDFA, self.sharedContextCache)
        self._predicates = None




    class Location_parseContext(ParserRuleContext):
        __slots__ = 'parser'

        def __init__(self, parser, parent:ParserRuleContext=None, invokingState:int=-1):
            super().__init__(parent, invokingState)
            self.parser = parser

        def location(self):
            return self.getTypedRuleContext(FeatureTableParser.LocationContext,0)


        def EOF(self):
            return self.getToken(FeatureTableParser.EOF, 0)

        def getRuleIndex(self):
            return FeatureTableParser.RULE_location_parse

        def enterRule(self, listener:ParseTreeListener):
            if hasattr( listener, "enterLocation_parse" ):
                listener.enterLocation_parse(self)

        def exitRule(self, listener:ParseTreeListener):
            if hasattr( listener, "exitLocation_parse" ):
                listener.exitLocation_parse(self)

        def accept(self, visitor:ParseTreeVisitor):
            if hasattr( visitor, "visitLocation_parse" ):
                return visitor.visitLocation_parse(self)
            else:
                return visitor.visitChildren(self)




    def location_parse(self):

        localctx = FeatureTableParser.Location_parseContext(self, self._ctx, self.state)
        self.enterRule(localctx, 0, self.RULE_location_parse)
        try:
            self.enterOuterAlt(localctx, 1)
            self.state = 22
            self.location()
            self.state = 23
            self.match(FeatureTableParser.EOF)
        except RecognitionException as re:
            localctx.exception = re
            self._errHandler.reportError(self, re)
            self._errHandler.recover(self, re)
        finally:
            self.exitRule()
        return localctx


    class Anticodon_parseContext(ParserRuleContext):
        __slots__ = 'parser'

        def __init__(self, parser, parent:ParserRuleContext=None, invokingState:int=-1):
            super().__init__(parent, invokingState)
            self.parser = parser

        def anticodon_value(self):
            return self.getTypedRuleContext(FeatureTableParser.Anticodon_valueContext,0)


        def EOF(self):
            return self.getToken(FeatureTableParser.EOF, 0)

        def getRuleIndex(self):
            return FeatureTableParser.RULE_anticodon_parse

        def enterRule(self, listener:ParseTreeListener):
            if hasattr( listener, "enterAnticodon_parse" ):
                listener.enterAnticodon_parse(self)

        def exitRule(self, listener:ParseTreeListener):
            if hasattr( listener, "exitAnticodon_parse" ):
                listener.exitAnticodon_parse(self)

        def accept(self, visitor:ParseTreeVisitor):
            if hasattr( visitor, "visitAnticodon_parse" ):
                return visitor.visitAnticodon_parse(self)
            else:
                return visitor.visitChildren(self)




    def anticodon_parse(self):

        localctx = FeatureTableParser.Anticodon_parseContext(self, self._ctx, self.state)
        self.enterRule(localctx, 2, self.RULE_anticodon_parse)
        try:
            self.enterOuterAlt(localctx, 1)
            self.state = 25
            self.anticodon_value()
            self.state = 26
            self.match(FeatureTableParser.EOF)
        except RecognitionException as re:
            localctx.exception = re
            self._errHandler.reportError(self, re)
            self._errHandler.recover(self, re)
        finally:
            self.exitRule()
        return localctx


    class LocationContext(ParserRuleContext):
        __slots__ = 'parser'

        def __init__(self, parser, parent:ParserRuleContext=None, invokingState:int=-1):
            super().__init__(parent, invokingState)
            self.parser = parser

        def location_term(self):
            return self.getTypedRuleContext(FeatureTableParser.Location_termContext,0)


        def accession(self):
            return self.getTypedRuleContext(FeatureTableParser.AccessionContext,0)


        def COLON(self):
            return self.getToken(FeatureTableParser.COLON, 0)

        def getRuleIndex(self):
            return FeatureTableParser.RULE_location

        def enterRule(self, listener:ParseTreeListener):
            if hasattr( listener, "enterLocation" ):
                listener.enterLocation(self)

        def exitRule(self, listener:ParseTreeListener):
            if hasattr( listener, "exitLocation" ):
                listener.exitLocation(self)

        def accept(self, visitor:ParseTreeVisitor):
            if hasattr( visitor, "visitLocation" ):
                return visitor.visitLocation(self)
            else:
                return visitor.visitChildren(self)




    def location(self):

        localctx = FeatureTableParser.LocationContext(self, self._ctx, self.state)
        self.enterRule(localctx, 4, self.RULE_location)
        self._la = 0 # Token type
        try:
            self.enterOuterAlt(localctx, 1)
            self.state = 31
            self._errHandler.sync(self)
            _la = self._input.LA(1)
            if _la==16:
                self.state = 28
                self.accession()
                self.state = 29
                self.match(FeatureTableParser.COLON)


            self.state = 33
            self.location_term()
        except RecognitionException as re:
            localctx.exception = re
            self._errHandler.reportError(self, re)
            self._errHandler.recover(self, re)
        finally:
            self.exitRule()
        return localctx


    class Location_termContext(ParserRuleContext):
        __slots__ = 'parser'

        def __init__(self, parser, parent:ParserRuleContext=None, invokingState:int=-1):
            super().__init__(parent, invokingState)
            self.parser = parser

        def COMPLEMENT(self):
            return self.getToken(FeatureTableParser.COMPLEMENT, 0)

        def LPAREN(self):
            return self.getToken(FeatureTableParser.LPAREN, 0)

        def location_term(self, i:int=None):
            if i is None:
                return self.getTypedRuleContexts(FeatureTableParser.Location_termContext)
            else:
                return self.getTypedRuleContext(FeatureTableParser.Location_termContext,i)


        def RPAREN(self):
            return self.getToken(FeatureTableParser.RPAREN, 0)

        def JOIN(self):
            return self.getToken(FeatureTableParser.JOIN, 0)

        def ORDER(self):
            return self.getToken(FeatureTableParser.ORDER, 0)

        def COMMA(self, i:int=None):
            if i is None:
                return self.getTokens(FeatureTableParser.COMMA)
            else:
                return self.getToken(FeatureTableParser.COMMA, i)

        def simple_location(self):
            return self.getTypedRuleContext(FeatureTableParser.Simple_locationContext,0)


        def getRuleIndex(self):
            return FeatureTableParser.RULE_location_term

        def enterRule(self, listener:ParseTreeListener):
            if hasattr( listener, "enterLocation_term" ):
                listener.enterLocation_term(self)

        def exitRule(self, listener:ParseTreeListener):
            if hasattr( listener, "exitLocation_term" ):
                listener.exitLocation_term(self)

        def accept(self, visitor:ParseTreeVisitor):
            if hasattr( visitor, "visitLocation_term" ):
                return visitor.visitLocation_term(self)
            else:
                return visitor.visitChildren(self)




    def location_term(self):

        localctx = FeatureTableParser.Location_termContext(self, self._ctx, self.state)
        self.enterRule(localctx, 6, self.RULE_location_term)
        self._la = 0 # Token type
        try:
            self.state = 53
            self._errHandler.sync(self)
            token = self._input.LA(1)
            if token in [2]:
                self.enterOuterAlt(localctx, 1)
                self.state = 35
                self.match(FeatureTableParser.COMPLEMENT)
                self.state = 36
                self.match(FeatureTableParser.LPAREN)
                self.state = 37
                self.location_term()
                self.state = 38
                self.match(FeatureTableParser.RPAREN)
                pass
            elif token in [3, 4]:
                self.enterOuterAlt(localctx, 2)
                self.state = 40
                _la = self._input.LA(1)
                if not(_la==3 or _la==4):
                    self._errHandler.recoverInline(self)
                else:
                    self._errHandler.reportMatch(self)
                    self.consume()
                self.state = 41
                self.match(FeatureTableParser.LPAREN)
                self.state = 42
                self.location_term()
                self.state = 47
                self._errHandler.sync(self)
                _la = self._input.LA(1)
                while _la==10:
                    self.state = 43
                    self.match(FeatureTableParser.COMMA)
                    self.state = 44
                    self.location_term()
                    self.state = 49
                    self._errHandler.sync(self)
                    _la = self._input.LA(1)

                self.state = 50
                self.match(FeatureTableParser.RPAREN)
                pass
            elif token in [14, 17]:
                self.enterOuterAlt(localctx, 3)
                self.state = 52
                self.simple_location()
                pass
            else:
                raise NoViableAltException(self)

        except RecognitionException as re:
            localctx.exception = re
            self._errHandler.reportError(self, re)
            self._errHandler.recover(self, re)
        finally:
            self.exitRule()
        return localctx


    class Simple_locationContext(ParserRuleContext):
        __slots__ = 'parser'

        def __init__(self, parser, parent:ParserRuleContext=None, invokingState:int=-1):
            super().__init__(parent, invokingState)
            self.parser = parser

        def NUMBER(self, i:int=None):
            if i is None:
                return self.getTokens(FeatureTableParser.NUMBER)
            else:
                return self.getToken(FeatureTableParser.NUMBER, i)

        def DOTDOT(self):
            return self.getToken(FeatureTableParser.DOTDOT, 0)

        def fuzzy_lt(self):
            return self.getTypedRuleContext(FeatureTableParser.Fuzzy_ltContext,0)


        def fuzzy_gt(self):
            return self.getTypedRuleContext(FeatureTableParser.Fuzzy_gtContext,0)


        def CARET(self):
            return self.getToken(FeatureTableParser.CARET, 0)

        def getRuleIndex(self):
            return FeatureTableParser.RULE_simple_location

        def enterRule(self, listener:ParseTreeListener):
            if hasattr( listener, "enterSimple_location" ):
                listener.enterSimple_location(self)

        def exitRule(self, listener:ParseTreeListener):
            if hasattr( listener, "exitSimple_location" ):
                listener.exitSimple_location(self)

        def accept(self, visitor:ParseTreeVisitor):
            if hasattr( visitor, "visitSimple_location" ):
                return visitor.visitSimple_location(self)
            else:
                return visitor.visitChildren(self)




    def simple_location(self):

        localctx = FeatureTableParser.Simple_locationContext(self, self._ctx, self.state)
        self.enterRule(localctx, 8, self.RULE_simple_location)
        self._la = 0 # Token type
        try:
            self.state = 68
            self._errHandler.sync(self)
            la_ = self._interp.adaptivePredict(self._input,5,self._ctx)
            if la_ == 1:
                self.enterOuterAlt(localctx, 1)
                self.state = 56
                self._errHandler.sync(self)
                _la = self._input.LA(1)
                if _la==14:
                    self.state = 55
                    self.fuzzy_lt()


                self.state = 58
                self.match(FeatureTableParser.NUMBER)
                self.state = 59
                self.match(FeatureTableParser.DOTDOT)
                self.state = 61
                self._errHandler.sync(self)
                _la = self._input.LA(1)
                if _la==15:
                    self.state = 60
                    self.fuzzy_gt()


                self.state = 63
                self.match(FeatureTableParser.NUMBER)
                pass

            elif la_ == 2:
                self.enterOuterAlt(localctx, 2)
                self.state = 64
                self.match(FeatureTableParser.NUMBER)
                self.state = 65
                self.match(FeatureTableParser.CARET)
                self.state = 66
                self.match(FeatureTableParser.NUMBER)
                pass

            elif la_ == 3:
                self.enterOuterAlt(localctx, 3)
                self.state = 67
                self.match(FeatureTableParser.NUMBER)
                pass


        except RecognitionException as re:
            localctx.exception = re
            self._errHandler.reportError(self, re)
            self._errHandler.recover(self, re)
        finally:
            self.exitRule()
        return localctx


    class AccessionContext(ParserRuleContext):
        __slots__ = 'parser'

        def __init__(self, parser, parent:ParserRuleContext=None, invokingState:int=-1):
            super().__init__(parent, invokingState)
            self.parser = parser

        def IDENTIFIER(self):
            return self.getToken(FeatureTableParser.IDENTIFIER, 0)

        def NUMBER(self):
            return self.getToken(FeatureTableParser.NUMBER, 0)

        def getRuleIndex(self):
            return FeatureTableParser.RULE_accession

        def enterRule(self, listener:ParseTreeListener):
            if hasattr( listener, "enterAccession" ):
                listener.enterAccession(self)

        def exitRule(self, listener:ParseTreeListener):
            if hasattr( listener, "exitAccession" ):
                listener.exitAccession(self)

        def accept(self, visitor:ParseTreeVisitor):
            if hasattr( visitor, "visitAccession" ):
                return visitor.visitAccession(self)
            else:
                return visitor.visitChildren(self)




    def accession(self):

        localctx = FeatureTableParser.AccessionContext(self, self._ctx, self.state)
        self.enterRule(localctx, 10, self.RULE_accession)
        self._la = 0 # Token type
        try:
            self.enterOuterAlt(localctx, 1)
            self.state = 70
            self.match(FeatureTableParser.IDENTIFIER)
            self.state = 73
            self._errHandler.sync(self)
            _la = self._input.LA(1)
            if _la==1:
                self.state = 71
                self.match(FeatureTableParser.T__0)
                self.state = 72
                self.match(FeatureTableParser.NUMBER)


        except RecognitionException as re:
            localctx.exception = re
            self._errHandler.reportError(self, re)
            self._errHandler.recover(self, re)
        finally:
            self.exitRule()
        return localctx


    class Fuzzy_ltContext(ParserRuleContext):
        __slots__ = 'parser'

        def __init__(self, parser, parent:ParserRuleContext=None, invokingState:int=-1):
            super().__init__(parent, invokingState)
            self.parser = parser

        def LT(self):
            return self.getToken(FeatureTableParser.LT, 0)

        def getRuleIndex(self):
            return FeatureTableParser.RULE_fuzzy_lt

        def enterRule(self, listener:ParseTreeListener):
            if hasattr( listener, "enterFuzzy_lt" ):
                listener.enterFuzzy_lt(self)

        def exitRule(self, listener:ParseTreeListener):
            if hasattr( listener, "exitFuzzy_lt" ):
                listener.exitFuzzy_lt(self)

        def accept(self, visitor:ParseTreeVisitor):
            if hasattr( visitor, "visitFuzzy_lt" ):
                return visitor.visitFuzzy_lt(self)
            else:
                return visitor.visitChildren(self)




    def fuzzy_lt(self):

        localctx = FeatureTableParser.Fuzzy_ltContext(self, self._ctx, self.state)
        self.enterRule(localctx, 12, self.RULE_fuzzy_lt)
        try:
            self.enterOuterAlt(localctx, 1)
            self.state = 75
            self.match(FeatureTableParser.LT)
        except RecognitionException as re:
            localctx.exception = re
            self._errHandler.reportError(self, re)
            self._errHandler.recover(self, re)
        finally:
            self.exitRule()
        return localctx


    class Fuzzy_gtContext(ParserRuleContext):
        __slots__ = 'parser'

        def __init__(self, parser, parent:ParserRuleContext=None, invokingState:int=-1):
            super().__init__(parent, invokingState)
            self.parser = parser

        def GT(self):
            return self.getToken(FeatureTableParser.GT, 0)

        def getRuleIndex(self):
            return FeatureTableParser.RULE_fuzzy_gt

        def enterRule(self, listener:ParseTreeListener):
            if hasattr( listener, "enterFuzzy_gt" ):
                listener.enterFuzzy_gt(self)

        def exitRule(self, listener:ParseTreeListener):
            if hasattr( listener, "exitFuzzy_gt" ):
                listener.exitFuzzy_gt(self)

        def accept(self, visitor:ParseTreeVisitor):
            if hasattr( visitor, "visitFuzzy_gt" ):
                return visitor.visitFuzzy_gt(self)
            else:
                return visitor.visitChildren(self)




    def fuzzy_gt(self):

        localctx = FeatureTableParser.Fuzzy_gtContext(self, self._ctx, self.state)
        self.enterRule(localctx, 14, self.RULE_fuzzy_gt)
        try:
            self.enterOuterAlt(localctx, 1)
            self.state = 77
            self.match(FeatureTableParser.GT)
        except RecognitionException as re:
            localctx.exception = re
            self._errHandler.reportError(self, re)
            self._errHandler.recover(self, re)
        finally:
            self.exitRule()
        return localctx


    class Anticodon_valueContext(ParserRuleContext):
        __slots__ = 'parser'

        def __init__(self, parser, parent:ParserRuleContext=None, invokingState:int=-1):
            super().__init__(parent, invokingState)
            self.parser = parser

        def LPAREN(self):
            return self.getToken(FeatureTableParser.LPAREN, 0)

        def POS(self):
            return self.getToken(FeatureTableParser.POS, 0)

        def COLON(self, i:int=None):
            if i is None:
                return self.getTokens(FeatureTableParser.COLON)
            else:
                return self.getToken(FeatureTableParser.COLON, i)

        def location(self):
            return self.getTypedRuleContext(FeatureTableParser.LocationContext,0)


        def COMMA(self, i:int=None):
            if i is None:
                return self.getTokens(FeatureTableParser.COMMA)
            else:
                return self.getToken(FeatureTableParser.COMMA, i)

        def AA(self):
            return self.getToken(FeatureTableParser.AA, 0)

        def amino_acid(self):
            return self.getTypedRuleContext(FeatureTableParser.Amino_acidContext,0)


        def SEQ(self):
            return self.getToken(FeatureTableParser.SEQ, 0)

        def sequence(self):
            return self.getTypedRuleContext(FeatureTableParser.SequenceContext,0)


        def RPAREN(self):
            return self.getToken(FeatureTableParser.RPAREN, 0)

        def getRuleIndex(self):
            return FeatureTableParser.RULE_anticodon_value

        def enterRule(self, listener:ParseTreeListener):
            if hasattr( listener, "enterAnticodon_value" ):
                listener.enterAnticodon_value(self)

        def exitRule(self, listener:ParseTreeListener):
            if hasattr( listener, "exitAnticodon_value" ):
                listener.exitAnticodon_value(self)

        def accept(self, visitor:ParseTreeVisitor):
            if hasattr( visitor, "visitAnticodon_value" ):
                return visitor.visitAnticodon_value(self)
            else:
                return visitor.visitChildren(self)




    def anticodon_value(self):

        localctx = FeatureTableParser.Anticodon_valueContext(self, self._ctx, self.state)
        self.enterRule(localctx, 16, self.RULE_anticodon_value)
        try:
            self.enterOuterAlt(localctx, 1)
            self.state = 79
            self.match(FeatureTableParser.LPAREN)
            self.state = 80
            self.match(FeatureTableParser.POS)
            self.state = 81
            self.match(FeatureTableParser.COLON)
            self.state = 82
            self.location()
            self.state = 83
            self.match(FeatureTableParser.COMMA)
            self.state = 84
            self.match(FeatureTableParser.AA)
            self.state = 85
            self.match(FeatureTableParser.COLON)
            self.state = 86
            self.amino_acid()
            self.state = 87
            self.match(FeatureTableParser.COMMA)
            self.state = 88
            self.match(FeatureTableParser.SEQ)
            self.state = 89
            self.match(FeatureTableParser.COLON)
            self.state = 90
            self.sequence()
            self.state = 91
            self.match(FeatureTableParser.RPAREN)
        except RecognitionException as re:
            localctx.exception = re
            self._errHandler.reportError(self, re)
            self._errHandler.recover(self, re)
        finally:
            self.exitRule()
        return localctx


    class Amino_acidContext(ParserRuleContext):
        __slots__ = 'parser'

        def __init__(self, parser, parent:ParserRuleContext=None, invokingState:int=-1):
            super().__init__(parent, invokingState)
            self.parser = parser

        def IDENTIFIER(self):
            return self.getToken(FeatureTableParser.IDENTIFIER, 0)

        def getRuleIndex(self):
            return FeatureTableParser.RULE_amino_acid

        def enterRule(self, listener:ParseTreeListener):
            if hasattr( listener, "enterAmino_acid" ):
                listener.enterAmino_acid(self)

        def exitRule(self, listener:ParseTreeListener):
            if hasattr( listener, "exitAmino_acid" ):
                listener.exitAmino_acid(self)

        def accept(self, visitor:ParseTreeVisitor):
            if hasattr( visitor, "visitAmino_acid" ):
                return visitor.visitAmino_acid(self)
            else:
                return visitor.visitChildren(self)




    def amino_acid(self):

        localctx = FeatureTableParser.Amino_acidContext(self, self._ctx, self.state)
        self.enterRule(localctx, 18, self.RULE_amino_acid)
        try:
            self.enterOuterAlt(localctx, 1)
            self.state = 93
            self.match(FeatureTableParser.IDENTIFIER)
        except RecognitionException as re:
            localctx.exception = re
            self._errHandler.reportError(self, re)
            self._errHandler.recover(self, re)
        finally:
            self.exitRule()
        return localctx


    class SequenceContext(ParserRuleContext):
        __slots__ = 'parser'

        def __init__(self, parser, parent:ParserRuleContext=None, invokingState:int=-1):
            super().__init__(parent, invokingState)
            self.parser = parser

        def IDENTIFIER(self):
            return self.getToken(FeatureTableParser.IDENTIFIER, 0)

        def getRuleIndex(self):
            return FeatureTableParser.RULE_sequence

        def enterRule(self, listener:ParseTreeListener):
            if hasattr( listener, "enterSequence" ):
                listener.enterSequence(self)

        def exitRule(self, listener:ParseTreeListener):
            if hasattr( listener, "exitSequence" ):
                listener.exitSequence(self)

        def accept(self, visitor:ParseTreeVisitor):
            if hasattr( visitor, "visitSequence" ):
                return visitor.visitSequence(self)
            else:
                return visitor.visitChildren(self)




    def sequence(self):

        localctx = FeatureTableParser.SequenceContext(self, self._ctx, self.state)
        self.enterRule(localctx, 20, self.RULE_sequence)
        try:
            self.enterOuterAlt(localctx, 1)
            self.state = 95
            self.match(FeatureTableParser.IDENTIFIER)
        except RecognitionException as re:
            localctx.exception = re
            self._errHandler.reportError(self, re)
            self._errHandler.recover(self, re)
        finally:
            self.exitRule()
        return localctx





