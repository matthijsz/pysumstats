import warnings


class SumStatsWarning(UserWarning):
    """Base class for SumStats warnings. Print behavior of these warnings can be adjusted using Pythons default warnings.simplefilter() method

    """
    pass


def sumstatswarn(message):
    warnings.warn(message, SumStatsWarning, 4)


# Force warnings.warn() to omit the source code line in the message for SumStatsWarnings
showwarning_orig = warnings.showwarning
warnings.showwarning = lambda message, category, filename, lineno, file=None, line=None: \
    showwarning_orig(message, category, filename, lineno, file=None, line='' if category == SumStatsWarning else None)
# Always print SumStatsWarning messages by default
warnings.simplefilter('always', SumStatsWarning)

