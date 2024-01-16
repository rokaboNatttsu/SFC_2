- [1. 目的意識と方針](#1-目的意識と方針)
- [2. モデルの仕様](#2-モデルの仕様)
  - [2.1. 定義式](#21-定義式)
  - [2.2. 行動方程式](#22-行動方程式)
    - [2.2.1. 価格決定](#221-価格決定)
    - [2.2.2. 解雇](#222-解雇)
    - [2.2.3. 雇用計画](#223-雇用計画)
    - [2.2.4. 生産、支出決定](#224-生産支出決定)
    - [2.2.5. 稼働率](#225-稼働率)
    - [2.2.6. 金利決定](#226-金利決定)
- [3. 明確なバグ](#3-明確なバグ)
- [4. 再現されている現象](#4-再現されている現象)

##  1. 目的意識と方針
いったん、複雑すぎになることを恐れず式をそろえてみる。それからモデルを完成させられそうなレベルまで単純化する

##  2. モデルの仕様
変数は2次元あるいは３次元の配列の形であらわす。構造体でエージェントを表すよりも良いと判断した

$X^e = \{\lambda_e X_{-1} + (1 - \lambda_e) X^e_{-1}\} \frac{X_{-1}}{X_{-2}}$ を使いたい。指数関数的な変化のトレンドに対して、直感的な値を出せると思う。ただ、期待値が線形写像ではなくなり、期待値同士の整合性は担保されなくなる。

$\sum_i \lambda_i = 1$ になるようなポートフォリオ配分をすべてのファンドと投資家と銀行が持つ。 $\lambda$ は適応的パラメータである

### 2.1. 定義式
- $btw(A, B, C)=min(max(A, B),C)$(ただし$A<C$)
- $YD_w = + w N_w + SS - T_{iw} - T_{ew} - i L_w$
- $YD_i = \Pi_{ci} + \Pi_{ki} + r_g GB_i - T_{ii} - T_{ei} - i_{bc} B_{ci} - i_{bk} B_{ki}$
- $N_w = N_c + N_k + N_b + N_f + N_g$
- $w = (w_c N_c + w_k N_k + w_b N_b + w_f N_f + w_g N_g)/N_w$
- $TC_c^e = pk^e I_c^e + w_c^e N_c^e + T_{vc}^e + T_{ec}^e + T_{cc}^e + i^e L_c - ibc^e B_c$
- $TC_k^e = (w_k^e N_k^e + T_{vk}^e + T_{ek}^e + T_{ck}^e + i^e L_k + i_{bk}^e B_k)$
- $UE=1$ if 失業している else $UE=0$(ABM限定)
- $UE=1 - \frac{N_w}{N}$(部門間モデル限定)
- $FN = \mu_{FN} + \sigma_{FN} * randn()$

### 2.2. 行動方程式
上から計算する
####    2.2.1. 価格決定
- $x \in \{c,k,b,f,g\}$
- $w_x = (1 - FN)w_{x-1}$ (if $\sum_{n=1}^{4} UE_{-n} > 2$ and $UE_{-1}=1$)(ABM限定)
- $w_x = (1 + FN)w_{x-1}$ (else if $\sum_{n=1}^{4} UE_{-n} <= 1$ and $UE_{-1}=0$)(ABM限定)
- $w_x = w_{x-1}$ (else)(ABM限定)
- $w_x = w_{x-1}\{1 + \lambda_{UE1}(\lambda_{UE2} - UE_{-1})\}$(部門間モデル限定)
- $mu_c = mu_{c-1}(1+FN)$ (if $\frac{IN}{S_{-1}} < \lambda_{mu}$)
- $mu_c = mu_{c-1}(1-FN)$ (else)
- $mu_k = mu_{k-1}(1+FN)$ (if $u_k^e > u^T$)
- $mu_k = mu_{k-1}(1-FN)$ (else)
- $p = (1 + mu_c)TC_c^e$
- $p_k = (1 + mu_k)TC_k^e$

####    2.2.2. 解雇
- $N_c-=\delta_c N_c; N_{UE} += 1$（部門間モデル限定）
- $N_k-=\delta_k N_k; N_{UE} += 1$（部門間モデル限定）
- $N_b-=\delta_b N_b; N_{UE} += 1$（部門間モデル限定）
- $N_f-=\delta_f N_f; N_{UE} += 1$（部門間モデル限定）
- $N_g-=\delta_g N_g; N_{UE} += 1$（部門間モデル限定）
- $UE = UE_{-1}$ (if $UE=1$)(ABM限定)
- $UE = 1$ (else if $rand() < \delta_x$)(ABM限定)
- $UE = 0$ (else)(ABM限定)

####    2.2.3. 雇用計画
- $N_c^D = btw((1-\beta_1)N_c, \gamma_c K_c, (1+\beta_2)N_c)$
- $N_k^D = btw((1-\beta_1)N_k, \gamma_k K_k, (1+\beta_2)N_k)$
- $N_b^D = btw((1-\beta_1)N_b, \gamma_b NL_b^e, (1+\beta_2)N_b)$
- $N_f^D = btw((1-\beta_1)N_f, \gamma_f NL_f^e, (1+\beta_2)N_f)$
- $N_g^D = btw((1-\beta_1)N_g, \gamma_g K_g, (1+\beta_2)N_g)$

- $N_{TMP} = \sum_{\{c,k,b,f,g\}} max(0, N_x^D-N_x)$
- while $N_{TMP} > 0$ and $N_{UE}>0$
  - $p_x=\frac{max(0, N_x^D-N_x)}{N_{TMP}}$
  - $r = rand()$
  - if $r < p_c$
    - $N_c += 1$
  - else if $r < p_c + p_k$
    - $N_k += 1$
  - else if $r < p_c + p_k + p_b$
    - $N_b += 1$
  - else if $r < p_c + p_k + p_b + p_f$
    - $N_f += 1$
  - else
    - $N_g += 1$
  - end
  - $N_{UE} -= 1$
  - $N_{TMP} -= 1$
- end

####    2.2.4. 生産、支出決定
- $S = btw\{0, C^e+C_g^e-IN+\eta_1 (C^e+C_g^e), min(\eta_2 N_c, \eta_3 K_c)\}$
- $I_k = -(\zeta_1 + \zeta_2 \frac{I_k}{K_k}) K_k + btw(0, (\frac{(\zeta_1 + \zeta_2 \frac{I_k}{K_k}) K_k + \{I_c^e + (\zeta_1 + \zeta_2 \frac{I_c^e}{K_c}) K_c\} + I_{gk}^e}{\eta_5 K_k} - u^T)K_k, min(\eta_4 N_k, \eta_5 K_k))$
- $p C_w = \alpha_1 YD_w^e + \alpha_2 NW_w$
- $p C_i = \alpha_3 YD_i^e + \alpha_4 NW_i$
- $C_g = \alpha_5 K_c$
- $I_c = -(\zeta_1 + \zeta_2 \frac{I_c}{K_c}) K_c + btw(0, (\frac{S}{\eta_3 K_c} - u^T)K_c, min(\eta_4 N_k, \eta_5 K_k, \eta_8 M_c + \eta_9 NL_c^e))$
- $I_g^D = \zeta_4 K_g$
- $I_{gk}^D = (1 - \zeta_3) \{I_g^D + (\zeta_1 + \zeta_2 \frac{I_g}{K_g}) K_g\}$
- $I_{gg}^D = \zeta_3 \{I_g^D + (\zeta_1 + \zeta_2 \frac{I_g}{K_g}) K_g\} - (\zeta_1 + \zeta_2 \frac{I_g}{K_g}) K_g$
- $I_{gkmax}=min(\eta_4 N_k, \eta_5 K_k) - [(\zeta_1 + \zeta_2 \frac{I_k}{K_k}) K_k + \{I_c + (\zeta_1 + \zeta_2 \frac{I_c}{K_c}) K_c\}]$
- $I_{gk} = btw(0, I_{gk}^D, min(I_{gkmax}, \eta_{10} M_k + \eta_{11} NL_k^e))$
- $I_{gg} = - (\zeta_1 + \zeta_2 \frac{I_g}{K_g}) K_g + btw(0, I_{gg}^D + (\zeta_1 + \zeta_2 \frac{I_g}{K_g}) K_g,\zeta_3 min(\eta_6 N_g, \eta_7 K_g) - (\zeta_1 + \zeta_2 \frac{I_g}{K_g} K_g))$
- $IN = IN_{-1} + S - C + Cg$

####  2.2.5. 稼働率
- $u_c = \frac{S}{\eta_3 K_c}$
- $u_k = \frac{I+(\zeta_1 + \zeta_2 \frac{I_k}{K_k}) K_k+I_{gk}}{\eta_5 K_k}$
- $u_g = \frac{I_{gg} + (\zeta_1 + \zeta_2 \frac{I_g}{K_g}) K_g}{\eta_7 K_g}$

####    2.2.6. 金利決定
ABMのバージョンは、売り注文と買い注文の数量と値段の決め方から作る必要がある。

##  3. 明確なバグ

##  4. 再現されている現象

