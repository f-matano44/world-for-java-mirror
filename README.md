# JA-WORLD
 | Download | Build status |
 |-----------|-----------|
 | [![Latest Release](https://gitlab.com/f-matano44/world-for-java/-/badges/release.svg)](https://gitlab.com/f-matano44/world-for-java/-/releases) | [![Build status](https://gitlab.com/f-matano44/world-for-java/badges/main/pipeline.svg)](https://gitlab.com/f-matano44/world-for-java/-/jobs) |


This is an independent Java port of WORLD vocoder (C++). <br>
[WORLD - a high-quality speech analysis, manipulation and synthesis system](https://github.com/mmorise/World) <br>

* JA-WORLD repository: https://gitlab.com/f-matano44/world-for-java
* mirroring: https://github.com/f-matano44/world-for-java


## Development environment and tools
* [Java 17 (LTS)](https://adoptium.net/temurin/releases/?version=17)
* [ant 1.10](https://ant.apache.org/bindownload.cgi)
* [VSCode](https://code.visualstudio.com/) + [Checkstyle for Java](https://marketplace.visualstudio.com/items?itemName=shengchen.vscode-checkstyle)
* [chatGPT (GPT-4.0)](https://chat.openai.com/)


## Usage
If you want to use some options, please read `test/TestApp.java`.

```java
import jp.f_matano44.ja_world.*;

double[] x = ... ;  // signal
int fs = ... ;      // sampling rate

double[][] f0_parameter = Dio.estimateF0(x, fs);
double[] _f0 = f0_parameter[0];
double[] t = f0_parameter[1];
double[] f0 = StoneMask.refineF0(x, _f0, t, fs);
double[][] sp = CheapTrick.estimateSp(x, f0, t, fs);
double[][] ap = D4C.estimateAp(x, f0, t, fs);

double[] y = Synthesis.getSignal(f0, sp, ap, fs);
```


## Calculation Error Between C++(original) and Java Implementations for `test/vaiueo2d.wav`

| Function | Maximum Difference |
|-----------|------------|
| Dio | 4.26e-13 |
| StoneMask | 5.68e-14 |
| CheapTrick | 4.81e-09 |
| D4C | 5.00e-10 |


## Build (to build/ directory)
```SH
ja-world/$ ant
```


## Run test application
```SH
ja-world/$ ant
ja-world/$ java -jar bin/TestApp.jar test/vaiueo2d.wav output.wav
```
