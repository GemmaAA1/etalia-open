from core.validators import CustomRegexValidator

def validate_feed_name(char):
    regex_validator = CustomRegexValidator(
        regex=r'^[\p{L}\s\']+$',
        message="Feed name must contain only letters",
        code='invalid')
    regex_validator(char)