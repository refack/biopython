from Bio.GenBank import AntlrParser

location_string = "join(1..100,200..300)"
location_obj = AntlrParser.parse_location(location_string)
print(location_obj)
assert location_obj == ('join', [('range', 1, 100, False, False), ('range', 200, 300, False, False)])

location_string = "complement(1..100)"
location_obj = AntlrParser.parse_location(location_string)
print(location_obj)
assert location_obj == ('complement', ('range', 1, 100, False, False))

print("All tests passed!")
