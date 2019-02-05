# liknorm

[![Travis](https://img.shields.io/travis/com/limix/liknorm.svg?style=flat-square&label=linux%20%2F%20macos%20build)](https://travis-ci.com/limix/liknorm) ![](https://img.shields.io/appveyor/ci/Horta/liknorm.svg?label=windows%20build&style=flat-square) [![Read the Docs (version)](https://img.shields.io/readthedocs/liknorm/latest.svg?style=flat-square)](http://liknorm.readthedocs.io/)

C library for computing moments of the product of an
exponential-family likelihood with a Normal distribution.

[Liknorm documentation](https://liknorm.readthedocs.io/).

## Install

You can install it via [conda](https://conda.io)

```bash
conda install -c conda-forge liknorm
```

Alternatively, one can compile and install it.
From Linux or MacOS systems, enter

```bash
bash <(curl -fsSL https://raw.githubusercontent.com/limix/liknorm/master/install)
```

from terminal. From Windows, enter

```dos
powershell -Command "(New-Object Net.WebClient).DownloadFile('https://raw.githubusercontent.com/limix/liknorm/master/install.bat', 'install.bat')" && install.bat
```

## Authors

- [Danilo Horta](https://github.com/horta)

## License

This project is licensed under the MIT License - see the [license file](https://raw.githubusercontent.com/limix/liknorm/master/LICENSE.md) for details.
