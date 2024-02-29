from setuptools import setup

__lib_name__ = "mactop3D"
__lib_version__ = "1.0.0"
__description__ = "deciphering chromatin domain, domain community and chromunity for 3D genome maps"
__url__ = "https://github.com/ydduanran/Mactop"
__author__ = "Ran Duan"
__author_email__ = "duanran@mail.ynu.edu.cn"
__license__ = "MIT"
__keywords__ = ["3D genome","Hi-C","Topologically associating domain (TAD)","TAD communities","Chromunities"]


setup(
    name = __lib_name__,
    version = __lib_version__,
    description = __description__,
    url = __url__,
    author = __author__,
    author_email = __author_email__,
    license = __license__,
    packages = ['mactop'],
    zip_safe = False
    )
