---
title: "Cointegration Analysis in Time Series"
author: "Kai Chien Chen"
format: 
  revealjs: 
    theme: "simple"
    transition: "slide"
    slide-number: true
    toc: true
    incremental: true
    keep-md: true
    math: true
    toc-depth: 1
execute: 
  echo: true

---




# 簡介

-   **共整合**（Cointegration）用於分析非平穩時間序列之間的長期關係。
-   常見於 **經濟學** 和 **金融學** 中。
-   範例：**購買力平價（PPP）** 與匯率。

## 時間序列與 I(1) 序列

-   **時間序列**：隨時間收集的數據點序列。
-   **I(1) 序列**：一階差分後成為定態的序列。
    -   例如：股票價格、GDP。

## I(1) 序列繪圖



::: {.cell}

```{.r .cell-code  code-fold="true"}
# 1. 生成一個非平穩的 I(1) 時間序列 (隨機漫步)
set.seed(123)  # 設置隨機種子，保證結果可重現
n <- 100
I1_series <- cumsum(rnorm(n))  # 累積隨機數生成 I(1) 序列

# 2. 對時間序列進行一階差分，應該變為定態
diff_series <- diff(I1_series)

# 3. 繪圖來比較原始 I(1) 序列和一階差分後的定態序列
plot(I1_series, type = 'l', col = 'blue', main = "I(1) 序列（非平穩）", ylab = "值")
lines(diff_series, type = 'l', col = 'red')
legend("topleft",lty = 1,col = c("blue","red"),legend = c("I(1)","一階差分後的 I(1)"))
```

::: {.cell-output-display}
![](midterm_presentation_files/figure-revealjs/unnamed-chunk-1-1.png){width=960}
:::
:::



## 範例

**購買力平價（PPP）** PPP 令 $q_t$ 為實際匯率的對數，$s_t$ 表示名義匯率的對數，$p_t$ 和 $p^*_{t}$ 分別表示國內與國外物價水平的對數。則有 $$q_t = s_t + p^*_{t}-p_t$$

# 什麼是共整合？

-   共整合：當兩個或多個 $I(1)$ 序列有一個**穩定的長期關係**時。
-   數學上，若存在一個向量 $A$ 使得：$A' y_t \sim I(0)$，其中 $y_t$ 是 I(1) 變量的向量。

# 共整合的限制

-   共整合假設 **長期穩定性**，但這在現實中不一定成立。
-   對於複雜數據，可能需要其他方法，例如 **非線性模型**。

# 誤差修正模型 {.smaller}

-   ECM 描述了共整合序列的**短期效果**和**長期效應**。

-   模型為：$\Delta y_t = \beta_0 \Delta Z_t + (\varphi - 1)[y_{t-1}-(\frac {\beta_0 + \beta_1}{1-\varphi})Z_{t-1}]+\epsilon_t$

-   其中：$$y_t = \varphi y_{t-1}+\beta_0 Z_t+\beta_1 Z_{t-1}+ \epsilon_t$$ $Z_t = Z_{t-1}+u_t$ 其中 $$u_t \sim N(0,\sigma^2)$$

# Engle-Granger 兩步驟法 (步驟 1)

1.  假設我們有兩個時間序列 $y_t$ 和 $z_t$。

2.  估計共整合：$$y_t = \beta_{0} + \beta_{1}z_{t} + e_t $$

3.  對 $\hat{e_t}$ 進行 ADF 檢驗：$$ \Delta \hat{e} = a_0 +a_{1}\hat{e_{t-1}} + \sum_{i=1}^{n}a_{i+1}\Delta \hat{e_{t-i}} + \epsilon_t $$

# Engle-Granger 兩步驟法 (步驟 2)

1. 零假設：
  $H_0:a_1=0$
  $H_1:a_1 < 0$
2. 如果拒絕 $H_0$，則兩者存在共整合。

# Johansen 檢驗

-   **Johansen 檢驗** 用於檢測多個共整合關係。

## VAR(p) 模型與 VECM 轉換

假設 $\mathbf{X}_t$ 為 $k$ 個變數的時間序列向量：
$$
\mathbf{X}_t = \mathbf{A}_1 \mathbf{X}_{t-1} + \mathbf{A}_2 \mathbf{X}_{t-2} + \dots + \mathbf{A}_p \mathbf{X}_{t-p} + \mathbf{\epsilon}_t
$$

將其改寫為差分形式的向量誤差修正模型（VECM）：
$$
\Delta \mathbf{X}_t = \Pi \mathbf{X}_{t-1} + \sum_{i=1}^{p-1} \Gamma_i \Delta \mathbf{X}_{t-i} + \mathbf{\epsilon}_t
$$

According to (Beveridge and Nelson JME 1981)

---

### 共整合矩陣 $\Pi$ 與長期關係

- **共整合矩陣** $\Pi = \alpha \beta^{'}$
  - $\beta$: 共整合向量矩陣，描述長期平穩關係
  - $\alpha$: 調整速度矩陣，表示變數回到均衡的速度

---

### $\Pi$ 的秩與共整合關係

根據 $\Pi$ 的秩，我們可以確定共整合向量的數量：

1. **秩為 0**（ $\text{rank}(\Pi) = 0$）：
   - 無共整合向量，即無長期均衡關係。
  
2. **秩介於 0 和 $k$ 之間**（$0 < \text{rank}(\Pi) < k$）：
   - 存在 $r$ 個共整合向量，代表 $r$ 個長期均衡關係。
  
3. **滿秩 $k$**（$\text{rank}(\Pi) = k$）：
   - 所有變數平穩，無需進行共整合檢定。

---

### $\Pi$的eigenvalue

1. 如果 
  $$rank(\Pi) = 0 \Rightarrow \lambda_{1} = \lambda_{2} = \cdots = \lambda_{k} = 0
  $$
  $$\Rightarrow log(1-\lambda_{i}) =0 \quad \forall i$$
  $$rank(\Pi) = k \Rightarrow log(1-\lambda_{i}) \not = 0 \quad \forall i
  $$
  則 $x_t$ 不存在共整合關係。

---

2. 如果 $rank(\Pi) = r$，且假設
$$
  \begin{cases}
\lambda_{1},\lambda_{2},\cdots,\lambda_{r} \not = 0 \\
\lambda_{r+1}=\lambda_{r+2}=\cdots=\lambda_{k}=0 \\
\end{cases}
$$
亦即，
$$
  \begin{cases}
log(1-\lambda_{i})\not = 0 \quad for \quad i=1,2,\cdots,r \\
log(1-\lambda_{i}) = 0 \quad for \quad i=r+1,r+2,\cdots,k \\
\end{cases}
$$
則 $x_{t}$ 存在共整合關係。

## Johansen 檢定步驟

1. 設定虛無假設與對立假設
   - $H_0$: 最大共整合階次為 $r$ (共整合關係數量)
   - $H_1$: 最大共整合階次為 $k$ or $r+1$
2. 計算兩種統計量：
   - **跡檢定統計量**：檢查 $r$ 是否足夠
   - **最大特徵值檢定統計量**：檢查 $r+1$ 個共整合向量的可能性

---

### 跡檢定統計量(Trace Test)

跡檢定的公式為：
$$
\text{Trace Statistic} = -T \sum_{i=r+1}^{k} \ln(1 - \lambda_i)
$$
- 若統計量 > 臨界值，拒絕 $H_0$，表示存在多於 $r$ 個共整合關係

---

### 最大特徵值檢定統計量(Max Test){.smaller}

最大特徵值檢定的公式為：
$$
\text{Max-Eigen Statistic} = -T \ln(1 - \lambda_{r+1})
$$

- 若統計量 > 臨界值，拒絕 $H_0$，進一步檢查 $r+1$ 個共整合向量。

- 較常使用。

- Step1: 檢定 $H_{0}:r=0 \quad vs. \quad H_{1}:r=1$ ，如果拒絕 $H_0 \rightarrow$ step2。

- Step2: 檢定 $H_{0}:r=1 \quad vs. \quad H_{1}:r=2$ ，如果拒絕 $H_0 \rightarrow$ $H_{0}:r=2 \quad vs. \quad H_{2}:r=3$，一直做到無法拒絕 $H_0$。

# Cointegration in R

-   **R package**: `urca` for conducting the Johansen test.

# R Code Example

## 1. Load Necessary Packages and data



::: {.cell}

```{.r .cell-code}
# 安裝必要packages
# install.packages("vars")
library(urca)
library(tidyverse)
library(tseries)
library(quantmod)
```
:::

::: {.cell result='hide'}

```{.r .cell-code}
# 下載資料
getSymbols("EWA", from="2006-04-26", to="2012-04-09")
```

::: {.cell-output .cell-output-stdout}

```
[1] "EWA"
```


:::

```{.r .cell-code}
getSymbols("EWC", from="2006-04-26", to="2012-04-09")
```

::: {.cell-output .cell-output-stdout}

```
[1] "EWC"
```


:::

```{.r .cell-code}
getSymbols("IGE", from="2006-04-26", to="2012-04-09")
```

::: {.cell-output .cell-output-stdout}

```
[1] "IGE"
```


:::

```{.r .cell-code}
# 提取調整後的價格
ewaAdj = unclass(EWA$EWA.Adjusted)
ewcAdj = unclass(EWC$EWC.Adjusted)
igeAdj = unclass(IGE$IGE.Adjusted)
```
:::



---

這三個資料為:

1. EWA - iShares MSCI Australia ETF
這是一個追蹤澳大利亞股票市場的 ETF，反映澳大利亞股票的整體表現。

2. EWC - iShares MSCI Canada ETF
這是一個追蹤加拿大股票市場的 ETF，反映加拿大市場的整體走勢。

3. IGE - iShares North American Natural Resources ETF
這是一個追蹤北美自然資源公司的 ETF，通常包括能源和材料類股，反映北美自然資源行業的表現。

## 2. Plot and describe the Data



::: {.cell layout-align="right"}
::: {.cell-output-display}
![](midterm_presentation_files/figure-revealjs/unnamed-chunk-4-1.png){fig-align='right' width=960}
:::
:::

::: {.cell layout-align="left"}
::: {.cell-output-display}


```{=html}
<div class="Rtable1"><table class="Rtable1">
<thead>
<tr>
<th class='rowlabel firstrow lastrow'></th>
<th class='firstrow lastrow'><span class='stratlabel'>Overall<br><span class='stratn'>(N=1499)</span></span></th>
</tr>
</thead>
<tbody>
<tr>
<td class='rowlabel firstrow'>EWA.Adjusted</td>
<td class='firstrow'></td>
</tr>
<tr>
<td class='rowlabel'>Mean (SD)</td>
<td>11.4 (2.28)</td>
</tr>
<tr>
<td class='rowlabel lastrow'>Median [Min, Max]</td>
<td class='lastrow'>12.1 [5.18, 15.7]</td>
</tr>
<tr>
<td class='rowlabel firstrow'>EWC.Adjusted</td>
<td class='firstrow'></td>
</tr>
<tr>
<td class='rowlabel'>Mean (SD)</td>
<td>19.5 (3.41)</td>
</tr>
<tr>
<td class='rowlabel lastrow'>Median [Min, Max]</td>
<td class='lastrow'>20.1 [9.99, 25.6]</td>
</tr>
<tr>
<td class='rowlabel firstrow'>IGE.Adjusted</td>
<td class='firstrow'></td>
</tr>
<tr>
<td class='rowlabel'>Mean (SD)</td>
<td>25.0 (4.57)</td>
</tr>
<tr>
<td class='rowlabel lastrow'>Median [Min, Max]</td>
<td class='lastrow'>24.8 [13.7, 35.5]</td>
</tr>
</tbody>
</table>
</div>
```


:::
:::




## 5. Augmented Dickey-Fuller Test (ADF)



::: {.cell result='hide'}

```{.r .cell-code}
# ADF檢定
adf.test(ewaAdj)
```

::: {.cell-output .cell-output-stdout}

```

	Augmented Dickey-Fuller Test

data:  ewaAdj
Dickey-Fuller = -1.732, Lag order = 11, p-value = 0.6918
alternative hypothesis: stationary
```


:::

```{.r .cell-code}
adf.test(ewcAdj)
```

::: {.cell-output .cell-output-stdout}

```

	Augmented Dickey-Fuller Test

data:  ewcAdj
Dickey-Fuller = -1.8042, Lag order = 11, p-value = 0.6612
alternative hypothesis: stationary
```


:::

```{.r .cell-code}
adf.test(igeAdj)
```

::: {.cell-output .cell-output-stdout}

```

	Augmented Dickey-Fuller Test

data:  igeAdj
Dickey-Fuller = -1.9476, Lag order = 11, p-value = 0.6005
alternative hypothesis: stationary
```


:::
:::



## 6. Select the Optimal Lag for VAR



::: {.cell}

```{.r .cell-code}
# 使用 VARselect 函數選擇最佳滯後階數
library(vars)
var_select <- VARselect(data.frame(ewaAdj, ewcAdj, igeAdj), lag.max = 10, type = "none")
var_select$selection
```

::: {.cell-output .cell-output-stdout}

```
AIC(n)  HQ(n)  SC(n) FPE(n) 
     4      2      2      4 
```


:::
:::



## 7. Johansen Cointegration Test



::: {.cell}

```{.r .cell-code  code-fold="false"}
# Johansen 共整合檢定
jotest.t <- ca.jo(data.frame(ewaAdj, ewcAdj, igeAdj), type="trace", K=3, ecdet="none", spec="longrun")
summary(jotest.t)
```

::: {.cell-output .cell-output-stdout}

```

###################### 
# Johansen-Procedure # 
###################### 

Test type: trace statistic , with linear trend 

Eigenvalues (lambda):
[1] 0.010707119 0.007278278 0.003526433

Values of teststatistic and critical values of test:

          test 10pct  5pct  1pct
r <= 2 |  5.28  6.50  8.18 11.65
r <= 1 | 16.21 15.66 17.95 23.52
r = 0  | 32.32 28.71 31.52 37.22

Eigenvectors, normalised to first column:
(These are the cointegration relations)

                EWA.Adjusted.l3 EWC.Adjusted.l3 IGE.Adjusted.l3
EWA.Adjusted.l3       1.0000000           1.000        1.000000
EWC.Adjusted.l3      -0.8294167        1913.905        5.868744
IGE.Adjusted.l3       0.1389803       -1583.616       -2.089788

Weights W:
(This is the loading matrix)

               EWA.Adjusted.l3 EWC.Adjusted.l3 IGE.Adjusted.l3
EWA.Adjusted.d    -0.001333851    7.233610e-06   -0.0007292815
EWC.Adjusted.d     0.030525548    7.509048e-06   -0.0009251048
IGE.Adjusted.d     0.042640500    1.751762e-05   -0.0006355543
```


:::
:::



---



::: {.cell}

```{.r .cell-code  code-fold="false"}
jotest.m <- ca.jo(data.frame(ewaAdj, ewcAdj, igeAdj), type="eigen", K=3, ecdet="none", spec="longrun")
summary(jotest.m)
```

::: {.cell-output .cell-output-stdout}

```

###################### 
# Johansen-Procedure # 
###################### 

Test type: maximal eigenvalue statistic (lambda max) , with linear trend 

Eigenvalues (lambda):
[1] 0.010707119 0.007278278 0.003526433

Values of teststatistic and critical values of test:

          test 10pct  5pct  1pct
r <= 2 |  5.28  6.50  8.18 11.65
r <= 1 | 10.93 12.91 14.90 19.19
r = 0  | 16.10 18.90 21.07 25.75

Eigenvectors, normalised to first column:
(These are the cointegration relations)

                EWA.Adjusted.l3 EWC.Adjusted.l3 IGE.Adjusted.l3
EWA.Adjusted.l3       1.0000000           1.000        1.000000
EWC.Adjusted.l3      -0.8294167        1913.905        5.868744
IGE.Adjusted.l3       0.1389803       -1583.616       -2.089788

Weights W:
(This is the loading matrix)

               EWA.Adjusted.l3 EWC.Adjusted.l3 IGE.Adjusted.l3
EWA.Adjusted.d    -0.001333851    7.233610e-06   -0.0007292815
EWC.Adjusted.d     0.030525548    7.509048e-06   -0.0009251048
IGE.Adjusted.d     0.042640500    1.751762e-05   -0.0006355543
```


:::
:::



## 8. Calculate and Plot Cointegrated Relationship



::: {.cell}

```{.r .cell-code  code-fold="false"}
# 提取共整合係數
alpha1 <- -0.8294165
alpha2 <-  0.1389800

# 計算共整合時間序列（CI）
CI <- coredata(ewaAdj) + alpha1 * coredata(ewcAdj) + alpha2 * coredata(igeAdj)

# 繪製共整合時間序列
plot(CI, type = "l", col = "purple", lwd = 2,
     main = "Cointegrated Relationship of EWA, EWC, and IGE",
     xlab = "Date", ylab = "Cointegrated Value")
abline(h = mean(CI), col = "darkgray", lty = 2) # 長期均衡水平
```

::: {.cell-output-display}
![](midterm_presentation_files/figure-revealjs/unnamed-chunk-10-1.png){width=960}
:::
:::



## 9. ADF Test for Cointegrated Series



::: {.cell}

```{.r .cell-code}
# ADF檢定共整合時間序列
adf.test(CI)
```

::: {.cell-output .cell-output-stdout}

```

	Augmented Dickey-Fuller Test

data:  CI
Dickey-Fuller = -3.1871, Lag order = 11, p-value = 0.09014
alternative hypothesis: stationary
```


:::

```{.r .cell-code}
kpss.test(CI)
```

::: {.cell-output .cell-output-stdout}

```

	KPSS Test for Level Stationarity

data:  CI
KPSS Level = 0.75332, Truncation lag parameter = 7, p-value = 0.01
```


:::
:::



## 10. Plot All Series and Cointegration



::: {.cell}

```{.r .cell-code  code-fold="true"}
# 繪製所有時間序列和共整合序列
plot(ewaAdj, type = "l", col = "darkgreen", lwd = 2, ylim = c(-3, 36),
     main = "Adjusted Closing Prices of ETFs (EWA, EWC, IGE)",
     xlab = "Date", ylab = "Adjusted Closing Price (USD)")

# 加入紅色線條表示 ewcAdj
lines(ewcAdj, col = "red", lwd = 2)

# 加入藍色線條表示 igeAdj
lines(igeAdj, col = "blue", lwd = 2)

# 加入紫色線條表示共整合序列
lines(CI, col = "purple", lwd = 2)

# 加入圖例
legend("topright", legend = c("EWA - Australia", "EWC - Canada", "IGE - North American Resources", "CI"),
       col = c("darkgreen", "red", "blue", "purple"), lwd = 2, cex = 0.8, box.lty = 0)
```

::: {.cell-output-display}
![](midterm_presentation_files/figure-revealjs/unnamed-chunk-12-1.png){width=960}
:::
:::



## 結果差異的原因

- trace 檢定統計量傾向於更保守，因為它是對所有特徵值的累積檢定，對所有的共整合向量的數量進行評估，檢查是否有額外共整合的可能性。
- 最大特徵值檢定更集中於每一個特徵值的增量效果，因此可能更敏感於檢測出一個特定的共整合關係。

## 為何結果不同

- 在結果中，trace 檢定檢測到較多的共整合向量數量（達到臨界值），而 max 檢定的顯著性水準較低，這可能說明：

  - 檢測力不同：trace 檢定在檢測多重共整合關係上較強，而 max 檢定更注重單一關係的顯著性。
  - 數據變動的影響：數據可能存在微弱的共整合關係，這在累積的 trace 檢定中更容易顯現，而在 max 檢定中可能無法通過嚴格的單一關係顯著性標準。

- 因此，trace 檢定和 max 檢定的不同結果表明數據中可能存在某些微弱的共整合結構，但其顯著性並非在每個檢定下都能確認。

# 結論{.smaller}


| **比較項目**                 | **Johansen 檢定**                             | **Engle-Granger 兩步驟法**                     |
|------------------------------|-----------------------------------------------|-----------------------------------------------|
| **目的**                     | 檢驗多變量系統中是否存在共整合關係                | 檢驗雙變量系統中的共整合關係                     |
| **模型形式**                 | 向量自回歸 (VAR) 模型                              | 兩步驟：第一步回歸檢驗單根，第二步檢驗誤差項的單根性   |
| **適用變數數量**             | 多變量（可檢驗多個變數間的共整合關係）              | 雙變量（通常僅限於檢驗兩個變數間的共整合）           |
| **優點**                     | 能同時檢驗多個共整合關係，適用於多變量系統             | 簡單易用，適合雙變量檢驗                              |
| **缺點**                     | 複雜且需要更大的樣本量，難度較高                       | 僅限於兩個變數，可能忽略多變量之間的共整合關係           |
| **適用條件**                 | 系統中的變數都需為同階整合                         | 變數必須是同階單根過程                                  |

---

| **比較項目**                 | **Trace Test**                                | **Max Test**                                  |
|------------------------------|-----------------------------------------------|-----------------------------------------------|
| **目標**                     | 檢驗在系統中是否存在至少 $r$ 個共整合向量         | 檢驗系統中是否存在精確 $r$ 個共整合向量            |
| **統計量**                   | $-T \sum_{i=r+1}^{k} \ln(1 - \lambda_i)$      | $-T \ln(1 - \lambda_{r+1})$                   |
| **解釋方式**                 | 檢驗所有特徵值的總和 | 主要關注最大特徵值，檢驗是否存在一個新的共整合關係     |
| **檢驗結果的靈敏度**         | 通常對多個共整合關係的檢測較為靈敏                 | 主要關注單一共整合關係的變化，較為精細                |

# 參考資料和markdown

- 林常青教授計量經濟(二)講義:Time Series Analysis (3): Spurious regression and Cointegration
- 時間序列分析: 總體經濟與財務金融之應用 / 陳旭昇 著
- https://www.quantstart.com/articles/Johansen-Test-for-Cointegrating-Time-Series-Analysis-in-R/
- [https://github.com/KaiChienChen/Statistical-Consulting/blob/main/test_presentation.qmd](https://github.com/KaiChienChen/Statistical-Consulting/blob/main/test_presentation.qmd)

# END
