# liknorm

[![Travis](https://img.shields.io/travis/limix/liknorm.svg?style=flat-square&label=linux%20%2F%20macos%20build)](https://travis-ci.org/limix/liknorm) [![AppVeyor](https://img.shields.io/appveyor/ci/Horta/liknorm.svg?style=flat-square&label=windows%20build)](https://ci.appveyor.com/project/Horta/liknorm)

C library for computing moments of the product of an
exponential-family likelihood with a Normal distribution.

## Install

You can install it via [conda](https://conda.io)

```bash
conda install -c conda-forge liknorm
```

Alternatively, one can compile and install it.
First, make sure you have [hcephes](https://github.com/limix/hcephes) installed.
Then, from Linux or MacOS systems, enter

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
