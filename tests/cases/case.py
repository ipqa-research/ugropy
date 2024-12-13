from typing import Union


# Class to define the test cases.
class Case:
    """Case class.

    This class is used to store the test cases for the ugropy models. The not
    given parameters (left as None) will not be checked.

    Parameters
    ----------
    identifier : str
        Identifier of the molecule.
    identifier_type : str
        Identifier type of the molecule (name, simles, mol, etc).
    cases_module : str
        Name of the module where the test case is stored. This is used only
        when the test case fails to show the dev where the test case is.
    commentary : str, optional
        Commentary of the test case. Usually to disctint something hard about
        the case. Also if a user found a bug could be mention here with the
        issue tag, by default "".
    unifac_result : Union[dict, None], optional
        Result of the UNIFAC model, by default None.
    psrk_result : Union[dict, None], optional
        Result of the PSRK model, by default None.
    joback_result : Union[dict, None], optional
        Result of the Joback model, by default None.
    """

    def __init__(
        self,
        identifier: str,
        identifier_type: str,
        cases_module: str,
        commentary: str = "",
        unifac_result: Union[dict, None] = None,
        psrk_result: Union[dict, None] = None,
        joback_result: Union[dict, None] = None,
    ):
        self.identifier = identifier
        self.identifier_type = identifier_type
        self.cases_module = cases_module
        self.commentary = commentary
        self.unifac_result = unifac_result
        self.psrk_result = psrk_result
        self.joback_result = joback_result
