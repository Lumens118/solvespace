This is the binary Software Development Kit (SDK) for Gmsh 4.15.0:

  * Operating system: Windows64-sdk (CYGWIN)
  * C++ compiler: /usr/bin/x86_64-w64-mingw32-g++.exe
  * C++ compiler ID: GNU
  * C++ compiler version: 11
  * C++ compiler flags: -fopenmp -O2 -g -DNDEBUG
  * Build options: 64Bit ALGLIB[contrib] ANN[contrib] Bamg Blas[petsc] Blossom Cgns DIntegration DomHex Eigen[contrib] Fltk Gmm[contrib] Hxt Jpeg Kbipack Lapack[petsc] MathEx[contrib] Med Mesh Metis[contrib] Mmg Mpeg Netgen Nii2mesh NoSocklenT ONELAB ONELABMetamodel OpenCASCADE OpenCASCADE-CAF OpenGL OpenMP OptHom PETSc Parser Plugins Png Post QuadMeshingTools QuadTri Solver TetGen/BR TinyXML2[contrib] Untangle Voro++[contrib] WinslowUntangler Zlib tinyobjloader

Gmsh is distributed under the terms of the GNU General Public License: see
share/doc/gmsh/LICENSE.txt and share/doc/gmsh/CREDITS.txt. For additional Gmsh
resources, see https://gmsh.info.

SDK layout:

  * lib/*gmsh*.{so,dll,dylib}: shared Gmsh library
  * lib/gmsh.dll.lib: import library (Windows only)
  * lib/gmsh.py: Python module
  * lib/gmsh.jl: Julia module
  * include/gmsh.h: C++ API header
  * include/gmshc.h: C API header
  * include/gmsh.h_cwrap: C++ wrapper of the C API (see the `Notes' below)
  * include/gmsh.f90: Fortran module
  * bin/gmsh: gmsh executable (linked with the shared Gmsh library)
  * share/doc/gmsh/tutorials/c++: C++ API tutorials
  * share/doc/gmsh/tutorials/c: C API tutorials
  * share/doc/gmsh/tutorials/python: Python API tutorials
  * share/doc/gmsh/tutorials/julia: Julia API tutorials
  * share/doc/gmsh/tutorials/fortran: Fortran API tutorials
  * share/doc/gmsh/examples/api: Other API examples

Notes:

  * The C API should work with most compilers.

  * The C++ API will only work if your compiler has an Application Binary
    Interface (ABI) that is compatible with the ABI of the compiler used to
    build this SDK (see above for the compiler ID and version).

  * If your C++ compiler does not have a compatible ABI and if there are no
    compatibility flags available, you can rename `gmsh.h_cwrap' as `gmsh.h':
    this implementation redefines the C++ API in terms of the C API. Using this
    header will lead to (slightly) reduced performance compared to using the
    native Gmsh C++ API from the original `gmsh.h' header, as it entails
    additional data copies between this C++ wrapper, the C API and the native
    C++ code.

    For example, the Windows SDK is currently compiled using the GNU Compiler
    Collection (GCC). To compile a C++ example with Microsoft Visual Studio 2017
    in the Visual Studio shell and run it, you would do:

    C:\gmsh-git-Windows64-sdk> ren include\gmsh.h gmsh.h_original
    C:\gmsh-git-Windows64-sdk> ren include\gmsh.h_cwrap gmsh.h
    C:\gmsh-git-Windows64-sdk> cl /Iinclude share\doc\gmsh\tutorials\c++\t1.cpp lib\gmsh.dll.lib
    C:\gmsh-git-Windows64-sdk> cd lib
    C:\gmsh-git-Windows64-sdk\lib> ..\t1.exe

  * To make it as portable and as easy-to-use as possible the shared Gmsh
    library embeds most Gmsh dependencies statically linked directly inside the
    share library: libgfortran, FreeType, OpenCASCADE, OpenBLAS, HDF5, MED,
    CGNS, FLTK, ... Linking your app with different versions of some of these
    dependencies can cause issues: in this case you should rebuild the Gmsh
    shared library from source and manage the dependencies consistently.

  * The shared Gmsh library also references shared system libraries. If some of
    these libraries are missing you will need to install them.


========================================================================================================================

这是 Gmsh 4.15.0 的二进制软件开发工具包（SDK）：

  * 操作系统：Windows64-sdk (CYGWIN)
  * C++ 编译器：/usr/bin/x86_64-w64-mingw32-g++.exe
  * C++ 编译器 ID：GNU
  * C++ 编译器版本：11
  * C++ 编译器标志：-fopenmp -O2 -g -DNDEBUG
  * 构建选项：64Bit ALGLIB[contrib] ANN[contrib] Bamg Blas[petsc] Blossom Cgns DIntegration DomHex Eigen[contrib] Fltk Gmm[contrib] Hxt Jpeg Kbipack Lapack[petsc] MathEx[contrib] Med Mesh Metis[contrib] Mmg Mpeg Netgen Nii2mesh NoSocklenT ONELAB ONELABMetamodel OpenCASCADE OpenCASCADE-CAF OpenGL OpenMP OptHom PETSc Parser Plugins Png Post QuadMeshingTools QuadTri Solver TetGen/BR TinyXML2[contrib] Untangle Voro++[contrib] WinslowUntangler Zlib tinyobjloader

Gmsh 根据 GNU 通用公共许可证条款分发：请参阅 share/doc/gmsh/LICENSE.txt 和 share/doc/gmsh/CREDITS.txt。
有关更多 Gmsh 资源，请参阅 https://gmsh.info。

SDK 布局：

  * lib/*gmsh*.{so,dll,dylib}：共享 Gmsh 库
  * lib/gmsh.dll.lib：导入库（仅限 Windows）
  * lib/gmsh.py：Python 模块
  * lib/gmsh.jl：Julia 模块
  * include/gmsh.h：C++ API 头文件
  * include/gmshc.h：C API 头文件
  * include/gmsh.h_cwrap：C API 的 C++ 包装器（见下面的"说明"）
  * include/gmsh.f90：Fortran 模块
  * bin/gmsh：gmsh 可执行文件（链接到共享 Gmsh 库）
  * share/doc/gmsh/tutorials/c++：C++ API 教程
  * share/doc/gmsh/tutorials/c：C API 教程
  * share/doc/gmsh/tutorials/python：Python API 教程
  * share/doc/gmsh/tutorials/julia：Julia API 教程
  * share/doc/gmsh/tutorials/fortran：Fortran API 教程
  * share/doc/gmsh/examples/api：其他 API 示例

说明：

  * C API 应该可以与大多数编译器一起工作。

  * C++ API 只有在您的编译器具有与用于构建此 SDK 的编译器的应用程序二进制接口（ABI）兼容的 ABI 时才能工作
    （请参阅上面的编译器 ID 和版本）。

  * 如果您的 C++ 编译器没有兼容的 ABI 并且没有可用的兼容性标志，您可以将 `gmsh.h_cwrap' 重命名为 `gmsh.h'：
    此实现根据 C API 重新定义 C++ API。与使用原始 `gmsh.h' 头文件中的原生 Gmsh C++ API 相比，使用此头文件
    会导致（略微）降低性能，因为它涉及此 C++ 包装器、C API 和原生 C++ 代码之间的额外数据复制。

    例如，Windows SDK 目前是使用 GNU 编译器集合（GCC）编译的。要在 Visual Studio shell 中使用 Microsoft
    Visual Studio 2017 编译 C++ 示例并运行它，您可以这样做：

    C:\gmsh-git-Windows64-sdk> ren include\gmsh.h gmsh.h_original
    C:\gmsh-git-Windows64-sdk> ren include\gmsh.h_cwrap gmsh.h
    C:\gmsh-git-Windows64-sdk> cl /Iinclude share\doc\gmsh\tutorials\c++\t1.cpp lib\gmsh.dll.lib
    C:\gmsh-git-Windows64-sdk> cd lib
    C:\gmsh-git-Windows64-sdk\lib> ..\t1.exe

  * 为了使其尽可能便携和易于使用，共享 Gmsh 库将大多数 Gmsh 依赖项静态链接直接嵌入到共享库中：
    libgfortran、FreeType、OpenCASCADE、OpenBLAS、HDF5、MED、CGNS、FLTK 等。将您的应用程序与这些依赖项的
    不同版本链接可能会导致问题：在这种情况下，您应该从源代码重新构建 Gmsh 共享库并一致地管理依赖项。

  * 共享 Gmsh 库还引用共享系统库。如果缺少某些库，您需要安装它们。
