- [1. 目的意識と方針](#1-目的意識と方針)
- [2. 定義式](#2-定義式)
- [3. 行動方程式](#3-行動方程式)
  - [3.1. 価格決定](#31-価格決定)
  - [3.2. 金利決定](#32-金利決定)
  - [3.3. 解雇](#33-解雇)
  - [3.4. 雇用計画](#34-雇用計画)
  - [3.5. 生産、支出決定](#35-生産支出決定)
  - [3.6. 稼働率](#36-稼働率)
  - [3.7. 租税(法人税と投資家の所得税以外)](#37-租税法人税と投資家の所得税以外)
  - [3.8. 利潤と法人税](#38-利潤と法人税)
  - [3.9. 投資家の所得税](#39-投資家の所得税)
  - [3.10. ポートフォリオ配分](#310-ポートフォリオ配分)
- [4. 恒等式](#4-恒等式)
  - [4.1. TFMの列の恒等式(モデル計算に使うやつ)](#41-tfmの列の恒等式モデル計算に使うやつ)
  - [4.2. TFMの行の恒等式(モデル計算に使うやつ)](#42-tfmの行の恒等式モデル計算に使うやつ)
  - [4.3. BSMの列の恒等式(モデル計算に使うやつ)](#43-bsmの列の恒等式モデル計算に使うやつ)
  - [4.4. BSMの行の恒等式(モデル計算に使うやつ)](#44-bsmの行の恒等式モデル計算に使うやつ)
  - [4.5. ストックとフローの接続の恒等式(モデル計算に使うやつ)](#45-ストックとフローの接続の恒等式モデル計算に使うやつ)
  - [4.6. 隠された恒等式](#46-隠された恒等式)
    - [4.6.1. 特に、初期値の整合性を担保するために使う恒等式](#461-特に初期値の整合性を担保するために使う恒等式)
- [5. 明確なバグ](#5-明確なバグ)
- [6. 再現されている現象](#6-再現されている現象)

##  1. 目的意識と方針
いったん、複雑すぎになることを恐れず式をそろえてみる。それからモデルを完成させられそうな水準まで単純化する。
- 金融市場から実物経済への影響の説明
- 内部資金と既存の資本の生産力に制限された、資本の再生産及び消費財の生産の性質の説明

を目指す

ABMの式に使う変数の添え字は正確に書いていない。書き直す必要がある

変数は2次元あるいは３次元の配列の形であらわす。構造体でエージェントを表すよりもコードが読み書きしやすくなると判断した

$X^e = \{\lambda_e X_{-1} + (1 - \lambda_e) X^e_{-1}\} \frac{X_{-1}}{X_{-2}}$ を使いたい。指数関数的な変化のトレンドに対して、直感的な値を出せると思う。ただ、期待値が線形写像ではなくなり、期待値同士の整合性は担保されなくなる。

$\sum_i r_i = 1$ になるようなポートフォリオ配分をすべてのファンドと投資家と銀行が持つ。 $r$ は適応的パラメータである

## 2. 定義式
- $X^e = \{\lambda_e X_{-1} + (1 - \lambda_e) X^e_{-1}\} \frac{X_{-1}}{X_{-2}}$ (ほかに明示しないときは、この方法で期待値を計算する)
- $btw(A, B, C)=min(max(A, B),C)$ (ただし $A < C$ )
- $YD_w = + w N_w + SS - T_{iw} - T_{ew} - i L_w$
- $YD_i = \Pi_{ci} + \Pi_{ki} + i_g GB_i - T_{ii} - T_{ei} - i_{bc} B_{ci} - i_{bk} B_{ki}$
- $N_w = N_c + N_k + N_b + N_f + N_g$
- $w = (w_c N_c + w_k N_k + w_b N_b + w_f N_f + w_g N_g)/N_w$
- $TC_c^e = pk^e I_c^e + w_c^e N_c^e + T_{vc}^e + T_{ec}^e + T_{cc}^e + i^e L_c - ibc^e B_c$
- $TC_k^e = (w_k^e N_k^e + T_{vk}^e + T_{ek}^e + T_{ck}^e + i^e L_k + i_{bk}^e B_k)$
- $UE=1$ if 失業している else $UE=0$(ABM限定)
- $UE=1 - \frac{N_w}{N}$(部門間モデル限定)
- $FN = \mu_{FN} + \sigma_{FN} * randn()$
- $Y = p(C+C_g) + p_k(I+I_{gg})$(部門間モデル限定)
- $Y = \sum pC_n + pC_g + \sum p_k (I_n + I_{ngk}) + p_k I_g - (\zeta_1 + \zeta_2 \frac{I_{gg}}{K_g} K_g)$(ABM限定)
- $E_c = E_{ci} + E_{cb} + E_{cf}$
- $E_k = E_{ki} + E_{kb} + E_{kf}$
- $NL_w = -p C_{ws} - (\zeta_1 + \zeta_2 \frac{C_{ws}}{K_w}) K_w + WS$
- $NL_i =  -p C_{is} - (\zeta_1 + \zeta_2 \frac{C_{is}}{K_i}) K_i + IS$
- $NL_c = -p_k I_c - p \Delta IN + \Pi_{cc}$
- $NL_k = -p_k I_k + \Pi_{kk}$
- $NL_b = w_b N_b - T_{cb} + \Pi_{cb} + \Pi_{kb} + i_g GB_b + i L + i_{bc} B_{cb} + i_{bk} B_{kb}$
- $NL_f = w_f N_f - T_{cf} + \Pi_{cf} + \Pi_{kf} + i_g GB_f - i L_f + i_{bc} B_{cf} + i_{bk} B_{kf}$
- $NL_g = -p_k I_{gg} + GS$
- $R_i = \frac{\Pi_{ci} + \Pi_{ki} + i_g GB_i + i_{bc} B_{ci} + i_{bk} B_{ki}}{M_i + p_{ci} E_{ci} + p_{ki} E_{ki} + B_{ci} + B_{ki}}$ 金融資本収益率
- $R_b = \frac{\Pi_{cb} + \Pi_{kb} + i_g GB_b + i L + i_{bc} B_{cb} + i_{bk} B_{kb}}{L + H_b + GB_b + p_{ec} E_{cb} + p_{ek} E_{kb} + B_{cb} + B_{kb}}$
- $R_f = \frac{\Pi_{cf} + \Pi_{kf} + i_g GB_f + i_{bc} B_{cf} + i_{bk} B_{kf}}{H_f + GB_f + p_{ec} E_{cf} + p_{ek} E_{kf} + B_{cf} + B_{kf}}$

## 3. 行動方程式
上から計算する
###    3.1. 価格決定
- $x \in \{ c,k,b,f,g \}$
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

###    3.2. 金利決定
ABMのバージョンは、売り注文と買い注文の数量と値段の決め方から作る必要がある。それから、計算順を検討したい
- $i = (1 - \iota_3)i_{-1} + \iota_3(\iota_1 + \iota_2 \frac{L}{Y^e})$(部門間モデル限定)
- $i_{bc} = (1 - \iota_3)i_{bc-1} + \iota_3(\iota_4 + \iota_5 \frac{B_c}{Y^e})$(部門間モデル限定)
- $i_{bk} = (1 - \iota_3)i_{bk-1} + \iota_3(\iota_4 + \iota_5 \frac{B_k}{Y^e})$(部門間モデル限定)
- $i_g$:外生的に定める
- $i = i_{-1}(1 + FN)$ (if $\frac{NL_b}{L} < \iota_6$ )(ABM限定)
- $i = i_{-1}(1 - FN)$ (else)(ABM限定)
- $i_{bc} = i_{bc-1}(1 + FN)$ (if $\frac{NL_c}{L_c + \Delta CP_c} < \iota_7$ )(ABM限定)
- $i_{bk} = i_{bc-1}(1 - FN)$ (else)(ABM限定)

###    3.3. 解雇
- $N_c-=\delta_c N_c; N_{UE} += 1$（部門間モデル限定）
- $N_k-=\delta_k N_k; N_{UE} += 1$（部門間モデル限定）
- $N_b-=\delta_b N_b; N_{UE} += 1$（部門間モデル限定）
- $N_f-=\delta_f N_f; N_{UE} += 1$（部門間モデル限定）
- $N_g-=\delta_g N_g; N_{UE} += 1$（部門間モデル限定）
- $UE = UE_{-1}$ (if $UE=1$)(ABM限定)
- $UE = 1$ (else if $rand() < \delta_x$)(ABM限定)
- $UE = 0$ (else)(ABM限定)

###    3.4. 雇用計画
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

###    3.5. 生産、支出決定
- $S = btw\{0, C^e+C_g^e-IN+\eta_1 (C^e+C_g^e), min(\eta_2 N_c, \eta_3 K_c)\}$
- $I_k = -(\zeta_1 + \zeta_2 \frac{I_k}{K_k}) K_k + btw(0, (\frac{(\zeta_1 + \zeta_2 \frac{I_k}{K_k}) K_k + \{I_c^e + (\zeta_1 + \zeta_2 \frac{I_c^e}{K_c}) K_c\} + I_{gk}^e}{\eta_5 K_k} - u^T)K_k, min(\eta_4 N_k, \eta_5 K_k))$
- $p C_{wf} = \lambda_{cw} (\alpha_1 YD_w^e + \alpha_2 NW_w)$
- $p C_{ws} = (1 - \lambda_{cw}) (\alpha_1 YD_w^e + \alpha_2 NW_w)$
- $p C_{if} = \lambda_{ci} (\alpha_3 YD_i^e + \alpha_4 NW_i)$
- $p C_{is} = (1 - \lambda_{cs}) (\alpha_3 YD_i^e + \alpha_4 NW_i)$
- $C_g = \alpha_5 K_c$
- $I_c = -(\zeta_1 + \zeta_2 \frac{I_c}{K_c}) K_c + btw(0, (\frac{S}{\eta_3 K_c} - u^T)K_c, min(\eta_4 N_k, \eta_5 K_k, \frac{\eta_8 M_c + \eta_9 NL_c^e}{p_k}))$
- $I_{gg}^D = \zeta_4 K_g$
- $I_{gk}^D = (1 - \zeta_3) \{I_{gg}^D + (\zeta_1 + \zeta_2 \frac{I_{gg}}{K_g}) K_g\}$
- $I_g^D = \zeta_3 \{I_{gg}^D + (\zeta_1 + \zeta_2 \frac{I_{gg}}{K_g}) K_g\} - (\zeta_1 + \zeta_2 \frac{I_{gg}}{K_g}) K_g$
- $I_{gkmax}=min(\eta_4 N_k, \eta_5 K_k) - [(\zeta_1 + \zeta_2 \frac{I_k}{K_k}) K_k + \{I_c + (\zeta_1 + \zeta_2 \frac{I_c}{K_c}) K_c\}]$
- $I_{gk} = btw(0, I_{gk}^D, min(I_{gkmax}, \frac{\eta_{10} M_k + \eta_{11} NL_k^e}{p_k}))$
- $I_g = - (\zeta_1 + \zeta_2 \frac{I_{gg}}{K_g}) K_g + btw(0, I_g^D + (\zeta_1 + \zeta_2 \frac{I_{gg}}{K_g}) K_g,\zeta_3 min(\eta_6 N_g, \eta_7 K_g) - (\zeta_1 + \zeta_2 \frac{I_{gg}}{K_g} K_g))$
- $IN = IN_{-1} + S - C + Cg$
- $SS = \lambda_{SS} \alpha_6 p C \frac{N_{UE}}{N} + (1-\lambda_{SS})SS_{-1}$

###  3.6. 稼働率
- $u_c = \frac{S}{\eta_3 K_c}$
- $u_k = \frac{I+(\zeta_1 + \zeta_2 \frac{I_k}{K_k}) K_k+I_{gk}}{\eta_5 K_k}$
- $u_g = \frac{I_g + (\zeta_1 + \zeta_2 \frac{I_{gg}}{K_g}) K_g}{\eta_7 K_g}$

###  3.7. 租税(法人税と投資家の所得税以外)
- $T_{vc} = \tau_1 p(C+C_g+\Delta IN)$
- $T_{vk} = \tau_2 p_k(I + I_{gk})$
- $T_{iw} = (\tau_3 + \tau_4 \frac{w N_w + SS}{Y}) (w N_w + SS)$(部門間モデル限定)
- $T_{iw} = (\tau_3 + \tau_4 \frac{(w + SS)N}{Y}) (w + SS)$(ABM限定)
- $T_{ew} = \tau_6 p K_w$
- $T_{ei} = \tau_6 p K_i$
- $T_{ec} = \tau_7 p_k K_c$
- $T_{ek} = \tau_7 p_k K_k$

###  3.8. 利潤と法人税
- $T_{cc} = max(0, \tau_8(p(C + C_g + \Delta IN) - w_c N_c - T_{vc} - T_{ec} - i L_c - i_{bc} B_c))$
- $\Pi_c = p(C + C_g + \Delta IN) - w_c N_c - T_{vc} - T_{ec} - T_{cc} - i L_c - i_{bc} B_c - T_{cc}$
- $\Pi_{cc} = (1 - \theta_c)\{\Pi_c - p_k I_c - p \Delta IN\}$
- $\Pi_{ci} = \theta_c\{\Pi_c - p_k I_c - p \Delta IN\}\frac{E_{ci}}{E_c}$
- $\Pi_{cb} = \theta_c\{\Pi_c - p_k I_c - p \Delta IN\}\frac{E_{cb}}{E_c}$
- $\Pi_{cf} = \theta_c\{\Pi_c - p_k I_c - p \Delta IN\}\frac{E_{cf}}{E_c}$
- $T_{ck} = max(0, \tau_8 (p_k(I + I_{gk}) - w_k N_k - T_{vk} - T_{ek} - i L_k - i_{bk} B_k))$
- $\Pi_k = p_k(I + I_{gk}) - w_k N_k - T_{vk} - T_{ek} - T_{ck} - i L_k - i_{bk} B_k - T_{ck}$
- $\Pi_{kk} = (1 - \theta_k)\{\Pi_k - p_k I_k\}$
- $\Pi_{ki} = \theta_k\{\Pi_k - p_k I_k\}\frac{E_{ki}}{E_k}$
- $\Pi_{kb} = \theta_k\{\Pi_k - p_k I_k\}\frac{E_{kb}}{E_k}$
- $\Pi_{kf} = \theta_k\{\Pi_k - p_k I_k\}\frac{E_{kf}}{E_k}$
- $T_{cb} = max(0, \tau_8(-w_b N_b - T_{cb} + \Pi_{cb} + \Pi_{kb} + i_g GB_b + i L + i_{bc} B_{cb} + i_{bk} B_{kb}))$
- $T_{cf} = max(0, \tau_8(-w_f N_f - T_{cf} + \Pi_{cf} + \Pi_{kf} + i_g GB_f - i L_f + i_{bc} B_{cf} + i_{bk} B_{kf}))$

###  3.9. 投資家の所得税
- $T_{ii} = \tau_5 (\Pi_{ci} + \Pi_{ki} + i_g GB_i + i_{bc} B_{ci} + i_{bk} B_{ki})$

###  3.10. ポートフォリオ配分
パラメータは $\kappa$ 、適応的に変化する値は $r$ 系列を使う
- $L_{w+1} = \kappa_1 YD_w^e$
- $E_{ci+1}=E_{c+1} \frac{NW_{i+1}-p K_{i+1}}{NW_{i+1}-p K_{i+1}+NW_{b+1}+NW_{f+1}}$(部門間モデル限定)
- $E_{cb+1}=E_{c+1} \frac{NW_{b+1}}{NW_{i+1}-p K_{i+1}+NW_{b+1}+NW_{f+1}}$(部門間モデル限定)
- $E_{cf+1}=E_{c+1} \frac{NW_{f+1}}{NW_{i+1}-p K_{i+1}+NW_{b+1}+NW_{f+1}}$(部門間モデル限定)
- $E_{ki+1}=E_{k+1} \frac{NW_{i+1}-p K_{i+1}}{NW_{i+1}-p K_{i+1}+NW_{b+1}+NW_{f+1}}$(部門間モデル限定)
- $E_{kb+1}=E_{k+1} \frac{NW_{b+1}}{NW_{i+1}-p K_{i+1}+NW_{b+1}+NW_{f+1}}$(部門間モデル限定)
- $E_{kf+1}=E_{k+1} \frac{NW_{f+1}}{NW_{i+1}-p K_{i+1}+NW_{b+1}+NW_{f+1}}$(部門間モデル限定)
- ミクロで、適応的な配分目標の更新をしたい。金融資産がだぶついたときの影響とか、金融不況や金融バブルの影響とか、軽全体の資産割合の変化とか、扱おうとすると、たぶんGodelyの行列とベクトルのやつは、静的すぎて、債権の目標保有額と存在する金額が一致しない問題を、別で解決する必要が出てくる。あるいは、債権価格が株式並みに価格変動する合理的な理由が必要になる
- $GB_{bi+1}^D = r_{i1} NW_{i+1}$
- $E_{ci+1}^D = \frac{r_{i2} NW_{i+1}}{p_{ec}}$
- $E_{ki+1}^D = \frac{r_{i3} NW_{i+1}}{p_{ek}}$
- $B_{ci+1}^D = r_{i4} NW_{i+1}$
- $B_{ki+1}^D = r_{i5} NW_{i+1}$
- ポートフォリオ目標/現在の保有残高/現在の株式市場価格から、売りオファーと買いオファーの金額と量を決める
  - 実際にアルゴリズムを書いてみて、良さそうな方法を採用したい。アルゴリズムのためしは、portfolio_simulation.ipynbで行う
- rの更新。 
  - $\sum_{m=1}^{M} (a^m R_m)$ ( $0 < a < 1, a^m$ の $m$ は指数)で収益率の高いエージェントの配分割合目標ベクトル $\vec{r_x}$ を、 $\vec{r_{+1}}=(1 - \lambda_r) \vec{r} + \lambda_r \vec{r_x}$ みたいにして中途半端にパクる
  - 預金を含めた資産配分割合の合計を１に維持しつつ、 $\vec{r}$ に摂動を加える

##   4. 恒等式
###  4.1. TFMの列の恒等式(モデル計算に使うやつ)
- $WS = -p C_{wf} + w N_w + SS - T_{iw} - T_{ew} - i L_w$
- $IS = -p C_{if} - T_{ii} - T_{ei} + \Pi_{ci} + \Pi_{ki} + i_g GB_i + i_{bc} B_{ci} + i_{bk} B_{ki}$
- $GS = -p C_g + p_k I_g - w_g N_g -SS + T_v + T_i + T_e + T_c - i_g GB$
- $\Delta M_w = NL_w + \Delta L_w$

###  4.2. TFMの行の恒等式(モデル計算に使うやつ)

###  4.3. BSMの列の恒等式(モデル計算に使うやつ)

###  4.4. BSMの行の恒等式(モデル計算に使うやつ)

###  4.5. ストックとフローの接続の恒等式(モデル計算に使うやつ)
- $\Delta L_w = L_{w+1} - L_w$
- $M_{w+1} = M_w + \Delta M_w$

###  4.6. 隠された恒等式
#### 4.6.1. 特に、初期値の整合性を担保するために使う恒等式

##  5. 明確なバグ

##  6. 再現されている現象

