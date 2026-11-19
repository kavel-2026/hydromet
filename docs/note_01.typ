#import "@local/modern-cug-report:0.1.3": *

#show: doc => template(doc, footer: "水文气象学 2026", header: "046241 段玉佯")

#align(center)[
  #text(size: 18pt, weight: "bold")[地转风公式的推导与 R 包实现]
  #v(0.2em)
  #text(size: 11pt)[班级：046241　姓名：段玉佯　学号：20241001286]
]

= 1 地转风的基本概念

== 1.1 定义

自由大气中，空气块在水平方向上受到水平气压梯度力与地转偏向力（科里奥利力）的作用。
当二者大小相等、方向相反时，空气做匀速直线运动，这种风称为*地转风*（geostrophic wind）。
地转风不穿越等压线，而是沿等压线吹，是天气图分析和大气环流研究的基本工具。

== 1.2 两种力

- *水平气压梯度力*：$bold(G) = -1/rho nabla p$，由高压指向低压，大小与气压梯度成正比；
- *地转偏向力*：$bold(C) = -f bold(k) times bold(V)$，垂直于风速，北半球指向运动方向的右侧，南半球指向左侧。

其中 $rho$ 为空气密度，$p$ 为气压，$bold(V) = (u, v)$ 为水平风速，$bold(k)$ 为垂直单位向量，
$f$ 为科里奥利参数。

== 1.3 地转平衡

当 $bold(G) + bold(C) = 0$ 时，水平气压梯度力与地转偏向力平衡，称为*地转平衡*。
此时空气块的水平加速度为零，风向、风速不再随时间变化。

= 2 地转平衡方程与求解

== 2.1 科里奥利参数

地转偏向力的大小与风速成正比，比例系数即科里奥利参数：

$ f = 2 Omega sin phi $ <eq_f>

其中 $Omega = 7.292115 times 10^(-5) "rad/s"$ 为地球自转角速度，$phi$ 为纬度。
北半球 $f > 0$，南半球 $f < 0$；赤道处 $f = 0$，地转平衡不成立。

== 2.2 水平运动方程

在 $f$ 平面、忽略摩擦力和水平加速度的条件下，水平运动方程退化为：

$ -f v = -1/rho (partial p)/(partial x), quad f u = -1/rho (partial p)/(partial y) $ <eq_balance>

式 #[@eq_balance] 表示水平气压梯度力与地转偏向力在 $x$、$y$ 两个方向上分别平衡。

== 2.3 地转风分量

由式 #[@eq_balance] 的第一式两边同除以 $-f$，得地转风的南北分量：

$ v_g = 1/(rho f) (partial p)/(partial x) $ <eq_v>

由第二式两边同除以 $f$，得地转风的东西分量：

$ u_g = -1/(rho f) (partial p)/(partial y) $ <eq_u>

这就是地转风的基本表达式：纬向风由经向气压梯度决定，经向风由纬向气压梯度决定，
风与气压梯度垂直。

== 2.4 矢量形式

把式 #[@eq_v]、#[@eq_u] 写成矢量形式：

$ bold(V)_g = 1/(rho f) bold(k) times nabla p $ <eq_vec>

叉乘 $bold(k) times nabla p$ 相当于把气压梯度矢量在北半球顺时针旋转 90°，
因此地转风沿等压线吹，且低压中心位于风的左侧。

= 3 公式的检验与物理意义

== 3.1 平衡关系的检验

把式 #[@eq_v]、#[@eq_u] 代回式 #[@eq_balance]，有

$ -f v_g = -1/rho (partial p)/(partial x), quad f u_g = -1/rho (partial p)/(partial y) $

两式恒成立；同时

$ bold(V)_g dot nabla p = u_g (partial p)/(partial x) + v_g (partial p)/(partial y) = 0 $

即地转风与气压梯度正交，这与“风沿等压线吹”的结论一致。

== 3.2 物理意义

- *风向*：地转风沿等压线吹。北半球背风而立，低压在左、高压在右；南半球相反。
- *风速*：$abs(bold(V)_g) = abs(nabla p) / (rho f)$。同样的气压梯度下，纬度越高 $f$ 越大，风速越小；
  赤道附近 $f -> 0$，风速趋于无穷，地转近似失效。
- *适用范围*：大尺度（水平尺度在千公里以上）、自由大气（摩擦层以上）、加速度较小的运动。
  近地面摩擦层、赤道与低纬、强曲率区域（如急流出口区）不宜使用地转近似。

= 4 R 包 hydromet 的实现

== 4.1 目录结构

R 包按“函数与脚本分离”的原则组织，模板来自
#link("https://github.com/rpkgs/rpkg.template")[rpkg.template]。本包的目录如下：

#figure(
  caption: [hydromet 包的目录结构],
  table(
    columns: (1fr, 2.2fr),
    align: left,
    table.header([路径], [内容]),
    [R/geostrophic.R], [地转风相关的 4 个函数],
    [R/Component.R], [去季节与分段回归函数（来自 RHtests）],
    [man/], [roxygen2 生成的帮助文档],
    [tests/testthat/], [testthat 单元测试],
    [docs/], [Typst 笔记与本作业],
    [.github/workflows/], [GitHub Actions 工作流],
  ),
)

== 4.2 函数设计

地转风部分的函数遵循“一次做好一件事”的原则，每个函数只完成一个计算：

#table(
  columns: (auto, 1fr),
  align: left,
  table.header([函数], [作用]),
  [`CoriolisParameter(lat)`], [计算 $f = 2 Omega sin phi$，输入纬度（度）],
  [`GeostrophicWind(dpdx, dpdy, f)`], [由气压梯度与 $f$ 计算 $(u_g, v_g)$],
  [`WindSpeed(u, v)`], [计算风速 $sqrt(u^2 + v^2)$],
  [`WindDirection(u, v)`], [计算气象风向（风的来向，0–360°）],
)

#block(breakable: false)[
  ```r
  GeostrophicWind <- function(dpdx, dpdy, f, rho = 1.2) {
    data.frame(
      u = -dpdy / (rho * f),
      v =  dpdx / (rho * f)
    )
  }
  ```
]

#pagebreak()

== 4.3 数值例子

取 45°N、空气密度 $rho = 1.2 "kg/m"^3$、纬向气压梯度 $10^(-3) "Pa/m"$（即 1 hPa/100 km）：

```r
f <- CoriolisParameter(45)          # 1.031e-04 s^-1
GeostrophicWind(dpdx = 1e-3, dpdy = 0, f = f)
```

```
  u        v
1 0 8.080724
```

即地转风为 8.08 m/s 的南风（$v > 0$，风向北吹，风向 180°）。
若改为经向气压梯度 $10^(-3) "Pa/m"$、纬向梯度为零，则 $u = -8.08$ m/s，
为西风（风向 270°）。可见北半球地转风相对气压梯度向右偏转 90°。

== 4.4 单元测试

测试用 testthat（第 3 版）编写，覆盖正常值、边界条件、错误输入和向量化：

#table(
  columns: (auto, auto, 1fr),
  align: left,
  table.header([测试文件], [用例数], [覆盖内容]),
  [`test-geostrophic.R`], [9], [科里奥利参数、地转平衡、南北半球、向量化、$f = 0$ 报错、风速风向],
  [`test-component.R`], [10], [去季节均值与异常、单/多断点回归、输入校验、SSE 接口],
)

在本地执行 `testthat::test_local()`，19 个用例、51 个断言全部通过：

```
component: .......................
geostrophic: ............................
DONE
```

== 4.5 持续集成

仓库配置了 4 个 GitHub Actions 工作流，推送后约 2–3 分钟即可看到结果：

#table(
  columns: (auto, auto, 1fr, auto),
  align: left,
  table.header([工作流], [运行环境], [内容], [结果]),
  [`R-CMD-check`], [windows-latest], [完整 `R CMD check`], [通过],
  [`R-CMD-check`], [ubuntu-latest], [完整 `R CMD check`], [通过],
  [`test-coverage`], [ubuntu-latest], [covr 覆盖率，上传 Codecov], [通过],
  [`pkgdown`], [ubuntu-latest], [构建并部署包文档网站], [未完成],
)

= 5 操作步骤

== 5.1 环境准备

1. R 4.6.2 已安装，但当前终端会话的 `PATH` 尚未刷新，`Rscript` 不可直接调用；
   用绝对路径 `C:\Program Files\R\R-4.6.2\bin\Rscript.exe` 完成全部命令。
2. 本机已有 roxygen2、knitr、magrittr，缺少 testthat 与 pkgload；改用清华 TUNA 镜像安装
   （官方源经代理下载需 27–44 s，镜像约 1 s）。
3. Typst 未安装：下载 `typst-x86_64-pc-windows-msvc.zip`（0.15.1），
   解压到 `D:\tools\typst`，并把该目录加入用户 `PATH`。
4. 模板依赖的字体 Times New Roman、Consolas、KaiTi、SimSun 均已安装，中文排版正常。

== 5.2 建立自己的仓库

课程要求“Use this template → Create a new repository”，而本地克隆的是老师的示例仓库
`kongdd/hello.R`，`origin` 没有写权限。因此新建 `kavel-2026/hydromet`，
修改 `DESCRIPTION` 中的包名、标题、作者，并把 `rpkg.template.Rproj` 重命名为 `hydromet.Rproj`。

== 5.3 添加函数与编写文档

1. 把 RHtests 的 `Component.R` 放入 `R/`：`RemoveSeasonal()`、`LS_OneBreak()`、
   `LS_MultiBreaks()` 与内部求解器 `lm_solve()`。
2. 新增 `R/geostrophic.R`，实现本文推导的 4 个函数。
3. 为每个函数补写 roxygen2 注释（`@param`、`@return`、`@examples`）；中文帮助文档需加
   `@encoding UTF-8`，否则 `checkRd()` 会报编码错误。
4. 执行 `roxygen2::roxygenise()` 自动生成 `man/` 与 `NAMESPACE`。

== 5.4 编写测试

在 `tests/testthat/` 下新建测试文件，用 `testthat::test_local()` 运行：
正常值校验解析解、边界条件（断点位置、$f = 0$）、错误输入（长度不一致）与向量化。
测试不比对硬编码的数字，而是校验物理关系：把地转风代回平衡方程、检验风与气压梯度正交、
南北半球符号相反；分段回归则用构造的解析解与 `lm()` 的结果作对照。

== 5.5 推送与查看结果

#block(breakable: false)[
  ```
  git add -A
  git commit -m "feat: 地转风函数、文档与测试"
  git push origin master
  ```
]

推送后 2–3 分钟即可在 Actions 页面看到结果。2026-11-19 22:53 的推送中，
`R-CMD-check` 用时 1.9 分钟（windows-latest 与 ubuntu-latest 均通过）、
`test-coverage` 用时 2.2 分钟并通过；`pkgdown` 需要现场编译 ragg、textshaping 等依赖，
本次运行 1 小时仍未完成，已取消，并把该工作流改为手动触发，避免每次推送都长时间占用 runner。
运行记录见
#link("https://github.com/kavel-2026/hydromet/actions")[仓库的 Actions 页面]。

= 6 存在的问题与处理

+ *R 不在当前会话的 PATH 中*：`Get-Command Rscript` 找不到命令，但用户环境变量里已有
  `C:\Program Files\R\R-4.6.2\bin`。原因是终端在修改 PATH 之前启动。
  处理：本作业用绝对路径调用；重启终端后恢复。

+ *R 包不完整*：只有 roxygen2、knitr、magrittr，缺少 testthat、pkgload。
  处理：从清华镜像安装；官方 CRAN 源经代理下载 `PACKAGES.gz` 需 27–44 s。

+ *未安装 Rtools*：`pkgbuild::has_rtools()` 返回 FALSE，无法编译含 C/C++ 的源码包。
  本包是纯 R 代码，`R CMD build`/`check` 不受影响；若以后加入 Rcpp 需补装 Rtools。

+ *Typst 未安装*：`typst` 命令不存在，模板 `@local/modern-cug-report:0.1.3` 无法编译。
  处理：下载官方 0.15.1 二进制，解压后加入 PATH。

+ *DESCRIPTION 被判为 Malformed*：Description 原以中文句号“。”结尾，R 不认为这是完整句子，
  `R CMD check` 报 “should contain one or more complete sentences”。
  处理：Description 改为英文整句。

+ *R CMD check 停在 checking package dependencies*：日志显示访问 bioconductor.org 的索引
  超时，前后约 10 分钟，看起来像卡死。处理：等待超时后检查继续；CI 网络正常，不受影响。

+ *克隆了老师的仓库*：`origin` 为 `kongdd/hello.R`，直接推送会 403。
  处理：新建自己的仓库并切换到 `hydromet` 包名。

+ *Git 直连 GitHub 超时*：`git push` 报 “Failed to connect to github.com:443 after 21068 ms”。
  本机需经系统代理 127.0.0.1:7897 才能访问。处理：`git -c http.proxy=... push`。

+ *没有 GitHub 凭据*：未安装 gh，Git Credential Manager 也没有缓存，推送时无法认证。
  处理：`git-credential-manager github login` 完成浏览器授权，再新建仓库并推送。

+ *中文 Rd 缺少编码声明*：`tools::checkRd()` 报大量
  “Non-ASCII contents without declared encoding”。
  处理：每个 roxygen 块加 `@encoding UTF-8`。

+ *覆盖率工作流会变红*：模板的 `test-coverage.yaml` 用 `secrets.CODECOV_TOKEN` 上传，
  学生仓库没有该 secret，推送时 `fail_ci_if_error` 为真会让工作流失败。
  处理：改为 `fail_ci_if_error: false`。

+ *pkgdown 部署缺少写权限*：GitHub 默认令牌只读，部署 `gh-pages` 会 403。
  处理：在工作流中加 `permissions: contents: write`。

+ *没有连接老师的服务器*：只在本地 Windows 与 GitHub Actions 的
  windows-latest、ubuntu-latest 上验证，未测试老师服务器上的 R 版本与包版本差异，
  这部分是本次作业的欠缺。

= 7 附录 A：Typst 使用自测

+ 标题用 `=` 的个数分级，不用 `#`。
+ 强调用 `*粗体*` 与 `_斜体_`；列表用 `-` 与 `+`。
+ 行内公式用 `$...$`，块公式单独成段并用 `$ ... $` 包裹。
+ 代码块用三个反引号包围，函数调用以 `#` 开头，如 `#figure(...)`。
+ 引入资源用 `#import "文件.typ"` 或 `#import "@包:版本": *`。

= 8 附录 B：编译与复现

```
cd docs
typst compile note_01.typ 046241_段玉佯.pdf
```

若需实时预览，可用 `typst compile --watch note_01.typ`，或在 VS Code 中使用 Tinymist 插件。
