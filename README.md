# MKPlugins

A collection of SuperCollider plugins by Mads Kjeldgaard.

Most of these are algorithms that I felt like exploring more in-depth at a C++ level to understand them better or simply found missing in the SuperCollider plugin landscape.


## Included plugins

See help files for more information on these, included links to sources/research where applicable:

**Chen**

Chen's chaotic double scroll attractor. Inspired by Bryan Head's implementation in the Mutable Instruments Stages alternative firmware.

**LPG**

A digital model of the Buchla Lowpass-Gate. Based on a paper/research by Julian Parker and Stefano D'Angelo.

### Requirements

- CMake >= 3.5
- SuperCollider source code

### Building

#### Preparations for a build
Clone the project:

    git clone https://github.com/madskjeldgaard/mkplugins
    cd mkplugins

This will get you the mkplugins.

Before continuing, you need two things: 

First, the [SuperCollider source code](github.com/supercollider/supercollider/) :

```bash
git clone github.com/supercollider/supercollider
```

Note where you placed this source code on your computer and copy the full path to it for use with CMake.

Secondly, you need the path to your SuperCollider extensions directory. This is where CMake will install the plugins.

You can get this by opening up SuperCollider and running this line of code

``` 
Platform.userExtensionDir
```

Note the path that it posts to your post window in SuperCollider, copy it and save it for the next part of the build process.


#### Actually Building
To build the plugins, run these commands. Note the two paths you need to put in to the CMake command, these are the ones you got from the preparatory steps above.

```bash
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE='Release' -DSC_PATH=/path/to/sc/source -DCMAKE_INSTALL_PREFIX=/path/to/extensions
cmake --build . --config Release
cmake --build . --config Release --target install
```


