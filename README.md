# JA-WORLD
This is WORLD that is remake to Java.<br>
https://github.com/mmorise/World<br>
JA-WORLD project is unrelated to original.<br>


## Development environment
* [Java 17 (LTS)](https://adoptium.net/temurin/releases/?version=17)
* [ant 1.10](https://ant.apache.org/bindownload.cgi)
* [VSCode](https://code.visualstudio.com/) + [Checkstyle for Java](https://marketplace.visualstudio.com/items?itemName=shengchen.vscode-checkstyle)
* [chatGPT (GPT-4.0)](https://chat.openai.com/)


## Usage
If you want to use some options, please read `test/TestApp.java`.

```java
import jp.f_matano44.ja_world.*;

double[][] f0Info = Dio.estimateF0(x, fs);
double[] _f0 = f0Info[0];
double[] t = f0Info[1];
double[] f0 = StoneMask.refineF0(x, _f0, t, fs);
double[][] sp = CheapTrick.estimateSp(x, f0, t, fs);
double[][] ap = D4C.estimateAp(x, f0, t, fs);

double[] y = Synthesis.getSignal(f0, sp, ap, fs);
```


## C++ vs Java: Maximum Difference in Computed Values for `test/vaiueo2d.wav`

| Function | Maximum Difference |
|-----------|------------|
| Dio | 3.41e-13 |
| StoneMask | 1.14e-13 |
| CheapTrick | 4.81e-09 |
| D4C | 5.00e-10 |


## Build (make jar file to build/ directory)
```SH
ja-world/$ ant
```


## Run test application
```SH
ja-world/$ ant
ja-world/$ java -jar bin/TestApp.jar test/vaiueo2d.wav output.wav
```
