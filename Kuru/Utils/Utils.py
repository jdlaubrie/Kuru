# COMPARE STRINGS WHICH MIGHT CONTAIN UNICODES
############################################################################
def insensitive(string):
    """Given a string, returns its lower/upper case insensitive string"""
    if getattr(str,'casefold',None) is not None:
        insen = lambda str_name: str_name.casefold()
    else:
        insen = lambda str_name: str_name.upper().lower()

    return insen(string)

