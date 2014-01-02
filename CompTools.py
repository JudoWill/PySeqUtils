__author__ = 'will'


def identity_score(a, b, weight=1, null=0):
    """Returns a particular value whenever it encounters identical inputs."""

    return weight if a.lower() == b.lower() else null

