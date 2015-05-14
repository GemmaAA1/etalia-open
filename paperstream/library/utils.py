from .constants import LANGUAGES, LANGUAGES_DETECT


def langcode_to_langpap(lang_code):
    """Convert the language code used in language detection in paperstream
    language code used in the app
    """
    return dict(LANGUAGES_DETECT)[lang_code.upper()]