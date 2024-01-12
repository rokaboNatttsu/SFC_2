##  目的意識と方針

##  モデルの仕様
変数は2次元あるいは３次元の配列の形であらわす。構造体でエージェントを表すよりも良いと判断した

$X^e = \{\lambda_e X_{-1} + (1 - \lambda_e) X^e_{-1}\} \frac{X_{-1}}{X_{-2}}$ を使いたい。指数関数的な変化のトレンドに対して、直感的な値を出せると思う。ただ、期待値が線形写像ではなくなり、期待値同士の整合性は担保されなくなる。

$\sum_i \lambda_i = 1$ になるようなポートフォリオ配分をすべてのファンドと投資家と銀行が持つ。 $\lambda$ は適応的パラメータである

### 定義式
$btw(A, B, C)=min(max(A, B),C)$(ただし$A<C$)
$YD_w = + w N_w + SS - T_{iw} - T_{ew} - i L_w$
$YD_i = \Pi_{ci} + \Pi_{ki} + r_g GB_i - T_{ii} - T_{ei} - i_{bc} B_{ci} - i_{bk} B_{ki}$
$N_w = N_c + N_k + N_b + N_f + N_g$
$w = (w_c N_c + w_k N_k + w_b N_b + w_f N_f + w_g N_g)/N_w$
$TC_c^e = pk^e I_c^e + w_c N_c^e + T_{vc}^e + T_{ec}^e + T_{cc}^e + i^e L_c - ibc^e B_c$
$UE=1$ if 失業している else $UE=0$(ABM限定)
$UE=1 - \frac{N_w}{N}$(部門間モデル限定)
$FN = \mu_{FN} + \sigma_{FN} * randn()$

### 行動方程式
上から計算する
####    価格決定
$w_c = (1 - FN)w_{c-1}$ (if $\sum_{n=1}^{4} UE_{-n} > 2$ and $UE_{-1}=1$)(ABM限定)
$w_c = (1 + FN)w_{c-1}$ (else if $\sum_{n=1}^{4} UE_{-n} <= 1$ and $UE_{-1}=0$)(ABM限定)
$w_c = w_{c-1}$ (else)(ABM限定)
$w_c = w_{c-1}\{1 + \lambda_{UE1}(\lambda_{UE2} - UE_{-1})\}$(部門間モデル限定)
$mu_c = mu_{c-1}(1+FN)$ (if $\frac{IN}{S_{-1}} < \lambda_{mu}$)
$mu_c = mu_{c-1}(1-FN)$ (else)
$p = (1 + mu_c)TC_c^e$

####    解雇
$N_c-=\delta_c N_c$
$N_k-=\delta_k N_k$
$N_b-=\delta_b N_b$
$N_f-=\delta_f N_f$
$N_g-=\delta_g N_g$

####    雇用計画
$N_c^D = btw((1-\beta_1)N_c, \gamma_c K_c, (1+\beta_2)N_c)$
$N_k^D = btw((1-\beta_1)N_k, \gamma_k K_k, (1+\beta_2)N_k)$
$N_b^D = btw((1-\beta_1)N_b, \gamma_b NL_b^e, (1+\beta_2)N_b)$
$N_f^D = btw((1-\beta_1)N_f, \gamma_f NL_f^e, (1+\beta_2)N_f)$
$N_g^D = btw((1-\beta_1)N_g, \gamma_g K_g, (1+\beta_2)N_g)$
完全雇用の枠を超えて求人が出るときは、給料が高いところから順に選択的に雇用が埋まることとする
$N_c = min(N_c^D, N_tmp); N_tmp -= N_c$
$N_k = min(N_k^D, N_tmp); N_tmp -= N_k$
$N_b = min(N_b^D, N_tmp); N_tmp -= N_b$
$N_f = min(N_f^D, N_tmp); N_tmp -= N_f$
$N_g = min(N_g^D, N_tmp); N_tmp -= N_g$

####    支出決定
$p C_w = \alpha_1 YD_w^e + \alpha_2 NW_w$
$p C_i = \alpha_3 YD_i^e + \alpha_4 NW_i$
$C_g = \alpha_5 K_c$

##  明確なバグ

##  再現されている現象

