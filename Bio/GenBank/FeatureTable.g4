grammar FeatureTable;

// Main entry points for parsing
location_parse: location EOF;
anticodon_parse: anticodon_value EOF;

// Location rules
location: (accession ':' )? location_term;

location_term
    : 'complement' '(' location_term ')'
    | ('join' | 'order') '(' location_term (',' location_term)* ')'
    | simple_location
    ;

simple_location
    : fuzzy_lt? NUMBER '..' fuzzy_gt? NUMBER
    | NUMBER '^' NUMBER
    | NUMBER
    ;

accession: IDENTIFIER ('.' NUMBER)?;
fuzzy_lt: '<';
fuzzy_gt: '>';

// Anticodon rules
anticodon_value: '(' 'pos' ':' location ',' 'aa' ':' amino_acid ',' 'seq' ':' sequence ')' ;

amino_acid: IDENTIFIER;
sequence: IDENTIFIER;

// Lexer rules
COMPLEMENT: 'complement';
JOIN: 'join';
ORDER: 'order';
POS: 'pos';
AA: 'aa';
SEQ: 'seq';

LPAREN: '(';
RPAREN: ')';
COMMA: ',';
COLON: ':';
DOTDOT: '..';
CARET: '^';
LT: '<';
GT: '>';

IDENTIFIER: [a-zA-Z_][a-zA-Z0-9_.]*;
NUMBER: [0-9]+;

WS: [ \t\r\n]+ -> skip;
