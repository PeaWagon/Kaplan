"""

Small version of PubChemPy.

"""

try:
    import requests
except ImportError:
    pass


PUG_REST_ERROR_CODES = {
    200: "Success",
    202: "Accepted (asynchronous operation pending)",
    400: "PUGREST.BadRequest: Request is improperly formed \
          (syntax error in the URL, POST body, etc.)",
    404: "PUGREST.NotFound: The input record was not found (e.g. invalid CID)",
    405: "PUGREST.NotAllowed: Request not allowed \
          (such as invalid MIME type in the HTTP Accept header)",
    504: "PUGREST.Timeout: The request timed out, from server overload or too broad a request",
    503: "PUGREST.ServerBusy: Too many requests or server is busy, retry later",
    501: "PUGREST.Unimplemented: The requested operation has \
          not (yet) been implemented by the server",
    500: ["PUGREST.ServerError: Some problem on the server side \
           (such as a database server down, etc.)",
          "PUGREST.Unknown: An unknown error occurred"],
}


# the space character is not included as this does
# not cause issues with PubChem
SPECIAL_CHARS = {
    "$", "&", "+", "=", ",", "/", ":", ";",
    "?", "@", "'", '"', "<", ">", "#", "%",
    "{", "}", "|", "\\", "^", "~", "[", "]",
    "`",
}


class WebError(Exception):
    """Exception class for web-related errors."""


def pubchem_request(input_value, input_type="name"):
    """Make HTTP request to PubChem database.

    Parameters
    ----------
    input_value : str
        Could be the name of the compound, cid, SMILES string,
        or InChI string etc.
    input_type : str
        Defaults to name. This is the input parameter type
        provided.

    Returns
    -------
    output as dictionary

    """
    valid_input_types = ["name", "inchikey", "inchi", "cid", "smiles"]
    assert input_type in valid_input_types

    # check internet is working
    try:
        requests.get("https://pubchemdocs.ncbi.nlm.nih.gov/")
    except NameError:
        raise WebError("Requests library is not installed.")
    except requests.exceptions.ConnectionError:
        raise WebError("Not connected to the internet.")

    # example URL:
    # https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/vioxx/property/InChI/TXT

    output_properties = {
        "CanonicalSMILES": None,
        "IsomericSMILES": None,
        "InChI": None,
        "InChIKey": None,
        "IUPACName": None,
        "Charge": None,
    }

    domain = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"

    # make sure input_value conforms with URL standards
    # otherwise, post request will be needed
    if any(c in str(input_value) for c in SPECIAL_CHARS):
        output_properties = pubchem_post(input_value, input_type, output_properties, domain)
    else:
        output_properties = pubchem_get(input_value, input_type, output_properties, domain)
    print(output_properties)
    try:
        output_properties["Charge"] = int(output_properties["Charge"])
    # TypeError: int() argument must be a string,
    # a bytes-like object or a number, not 'NoneType'
    except TypeError:
        pass
    return output_properties


def pubchem_get(input_value, input_type, output_properties, domain):
    """Use HTTP GET request to query PubChem."""
    for property in output_properties:
        url = f"{domain}/{input_type}/{input_value}/property/{property}/TXT"
        result = requests.get(url)
        if result.status_code == 200:
            content = result._content.decode().split("\n")
            output_properties[property] = content[0]

    cid = requests.get(f"{domain}/{input_type}/{input_value}/cids/TXT")
    if cid.status_code == 200:
        content = cid._content.decode().split("\n")
        output_properties["cid"] = content[0]
    else:
        output_properties["cid"] = None

    sdf_3d = requests.get(f"{domain}/{input_type}/{input_value}/SDF?record_type=3d")
    if sdf_3d.status_code == 200:
        output_properties["sdf"] = sdf_3d._content.decode()
    else:
        output_properties["sdf"] = None
    return output_properties


def pubchem_post(input_value, input_type, output_properties, domain):
    """Use HTTP POST request to query PubChem.

    Notes
    -----
    Required for pubchem requests using InChI or
    SMILES/name with special characters.

    """

    header = {"Content-Type": "application/x-www-form-urlencoded"}
    data_input = {input_type: str(input_value)}

    for property in output_properties:
        url = f"{domain}/{input_type}/property/{property}/TXT"
        result = requests.post(url, headers=header, data=data_input)
        if result.status_code == 200:
            # strip newline character
            content = result._content.decode().split("\n")
            output_properties[property] = content[0]

    cid = requests.post(
        f"{domain}/{input_type}/cids/TXT",
        headers=header, data=data_input
    )
    if cid.status_code == 200:
        content = cid._content.decode().split("\n")
        output_properties["cid"] = content[0]
    else:
        output_properties["cid"] = None

    sdf_3d = requests.post(
        f"{domain}/{input_type}/SDF?record_type=3d",
        headers=header, data=data_input
    )
    if sdf_3d.status_code == 200:
        output_properties["sdf"] = sdf_3d._content.decode()
    else:
        output_properties["sdf"] = None
    return output_properties
