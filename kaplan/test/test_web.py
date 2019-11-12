"""

Tests for the web module from Kaplan.

"""

from kaplan.web import pubchem_request


def test_pubchem_request():
    """Test pubchem_request function from web module."""
    # how many results should be in each output dictionary
    key_count = 8

    # test cases that should all return values
    tests = ["butane", "ethyl alcohol", "benzene", "caffeine", "HCN"]
    for t in tests:
        result = pubchem_request(t)
        assert len(result) == key_count
        for value in result.values():
            assert value is not None

    result = pubchem_request(2244, "cid")
    assert len(result) == key_count
    for value in result.values():
        assert value is not None

    result = pubchem_request("CCCC", "smiles")
    assert len(result) == key_count
    for value in result.values():
        assert value is not None

    result = pubchem_request("ISAKRJDGNUQOIC-UHFFFAOYSA-N", "inchikey")
    assert len(result) == key_count
    for value in result.values():
        assert value is not None

    result = pubchem_request(
        ("InChI=1S/C46H70O2/c1-10"
         "-11-12-20-37(2)21-13-22-38(3)23-14-24-39(4)25"
         "-15-26-40(5)27-16-28-41(6)29-17-30-42(7)31-18"
         "-32-43(8)35-36-44-33-19-34-45(48-9)46(44)47/h"
         "10-11,19,21,23,25,27,29,31,33-35,47H,12-18,20"
         ",22,24,26,28,30,32,36H2,1-9H3"), "inchi"
    )
    assert len(result) == key_count
    for key, value in result.items():
        if key != "sdf":
            assert value is not None

    result = pubchem_request(
        "InChI=1S/C10H16O/c1-7(2)9-5-4-8(3)6-10(9)11/h8-9H,1,4-6H2,2-3H3",
        "inchi"
    )
    assert len(result) == key_count
    for value in result.values():
        assert value is not None

    result = pubchem_request("invalid molecule")
    assert len(result) == key_count
    for value in result.values():
        assert value is None

    # a 3D record for this compound does not exist; therefore, get 404
    # not found error with 3D request
    # 2D request should work
    result = pubchem_request("Tetracontane")
    assert len(result) == key_count
    for key, value in result.items():
        if key != "sdf":
            assert value is not None

    no3d = []
    for i in range(100):
        result = pubchem_request(i, "cid")
        if not result["sdf"]:
            no3d.append(i)
