# SPEC.md — genesis-spacefix の変更点

このリポジトリは GENESIS の `spdyn` に、原子を空間に固定する機能(space fix / 固定原子)を追加したものです。
このファイルは、その変更が何をしていて、どこに手を入れたかをまとめたものです。

関連コミット(`dev` ブランチ):

- `9fbbf5a` Add fixed atoms (fixatm_select_index) for spdyn REMD (gREST)
- `e0bdda2` Accept fixed atoms by atom number (fixatm_atom_index)
- `7b97805` Allow SHAKE/SETTLE with fixed atoms in spdyn REMD (AMBER)

## 1. 何をする機能か

spdyn の REMD(gREST を含む)実行で、指定した原子を空間に固定したまま MD を走らせる機能。
固定した原子は

- 質量・逆質量が 0 になり、力を受けても速度が変わらない。
- 初期速度も 0 にする。
- 自由度(DOF)の計算から除外する。
- 重心(COM)の除去の対象から外す。

という形で「動かない」原子として扱われる。

atdyn にはもともと `fixatm_select_index`(`[CONSTRAINTS]`)という同種の機能があるが、spdyn には無かった。
今回の実装は、その考え方を spdyn のドメイン分割・MPI 並列の構造に合わせて作り直したもの。

## 2. 対応範囲

| 項目 | 対応状況 |
|---|---|
| REMD(gREST を含む、`setup_spdyn_remd`) | 対応 |
| MD / 最小化 / RPATH | 未対応(固定原子指定があるとエラーで停止) |
| 並列 I/O(PIO、再開ファイル名が `rst(...)` 形式のもの) | 未対応(エラーで停止) |
| FEP(`[ALCHEMY]`) | 未対応(エラーで停止) |
| 水分子を固定原子に含める | 不可(エラーで停止、§5 参照) |
| `rigid_bond = YES`(SHAKE/SETTLE)との併用 | **AMBER のみ対応**。CHARMM・GROMACS はエラーで停止 |

## 3. 入力キーワード(`[CONSTRAINTS]`)

```
[CONSTRAINTS]
fixatm_atom_index   = 1-94 200 350-360   # 原子番号(1始まり)。a-b の範囲指定可、空白区切り
fixatm_select_index = 1 2                # [SELECTION] のグループ番号(複数可)
```

- 両方指定した場合は和集合(重複は 1 回だけ数える)。
- どちらも空なら、固定原子機能そのものが無効になる(通常の spdyn と同じ挙動)。
- 水分子(残基名が `TIP`/`WAT`/`MIZ`/`SOL` で始まるもの)を選択するとエラーで停止する。水は専用の積分経路(テーブル質量)を使うため、質量を 0 にしても効かないためこの検査を入れている。

## 4. 実装の仕組み

### 4.1 固定原子の選択(`setup_fixatm_spdyn`、`src/spdyn/sp_constraints.fpp`)

- `fixatm_atom_index` をパース(自前の空白分割。`a-b` の範囲は自前で展開する。`split_string` は lib に実体が無くリンクできなかったため使っていない)。
- `fixatm_select_index` のグループ番号を `select_atom`(lib の選択式パーサ)で原子に変換。
- 水分子を含んでいればここでエラー。
- 戻り値は、分子全体サイズの `logical` 配列 `fixed(:)`(どの原子が固定か)と、固定原子数 `num_fixatm`。

呼び出し元: `setup_spdyn_remd`(`src/spdyn/sp_setup_spdyn.fpp`)。`setup_domain` より前に呼ぶ。

### 4.2 質量・速度・自由度(`apply_fixatm_domain`、`sp_constraints.fpp` / 呼び出しは `sp_setup_spdyn.fpp`)

質量を `molecule` レベルで 0 にすると、水の判定(`Setup_Solute_And_Water> # of water is incorrect`)が壊れることが分かったため、**`setup_domain` より後**、ドメイン分割された `domain` 配列に対して適用する。

```
domain%mass(ix,i)     = 0
domain%inv_mass(ix,i) = 0
domain%velocity(1:3,ix,i) = 0
```

`inv_mass` は原子が別セルへ移動するたびに `mass` から再計算される(`sp_communicate.fpp`、`sp_domain.fpp` の複数箇所)。既存コードは `if (mass(ix,i) > EPS) inv_mass = 1/mass` という形で書かれているため、`mass = 0` を一度設定しておけば、再分配後もガードに引っかかって `inv_mass = 0` のままになる。これが「固定」の実体。

自由度は、`setup_domain` が重複数(`domain%num_duplicate`、通常は 1)を掛けた後に

```
update_num_deg_freedom('After setup of fixed atoms', -3*num_fixatm*domain%num_duplicate, molecule%num_deg_freedom)
```

で減らす。

### 4.3 重心(COM)除去からの除外(`src/spdyn/sp_assign_velocity.fpp`)

`stop_trans_rotation` の並進・回転除去の各ループに

```
if (mass(ix,i) <= EPS) cycle
```

を追加。質量 0 の原子(=固定原子)を COM の計算・速度補正の対象から外す。これを入れないと、COM 除去で固定原子の速度が変わってしまい、座標が動いてしまう。

### 4.4 SHAKE/SETTLE との併用(AMBER のみ)

固定原子を含む水素結合グループ(重原子 + その水素)が SHAKE の対象になっていると、固定原子と自由原子の距離を一定に保とうとして自由原子側が引っ張られる、あるいは固定原子自体が動かされる。これを避けるため、固定原子を含むグループの水素は「軽い原子(light atom)」のフラグを外し、SHAKE の対象から外す。

- `unconstrain_fixed_hydrogen`(`src/spdyn/sp_domain.fpp`): `setup_domain` の `check_light_atom_name` の直後に呼ぶ。分子の全結合を見て、「水素原子とその重原子パートナーのどちらかが固定原子」であるグループを特定し、そのグループに属する水素の `molecule%light_atom_name` / `molecule%light_atom_mass` を `.false.` にする。
  - 軽い原子でなくなった水素は、以後の処理(H グループの作成、ドメイン振り分け、通信、移動)で「ふつうの重原子」と同じに扱われる。これにより、グループごと SHAKE の対象から外れ、かつドメイン振り分けからも漏れない(漏れると原子が計算から消えるという不具合になることを実機で確認済み。§7 参照)。
- `src/spdyn/sp_enefunc_amber.fpp`(`setup_enefunc_bond_constraint`): 上記で「軽い原子」でなくなった水素の結合は、SHAKE グループ(`HGr`)に一致しなくなる。従来はそのような H-bond は無視されて結合項が存在しないままだったため、`found_hgr` フラグを追加し、一致しなかった場合はその結合を**通常の結合項**として `enefunc%bond_list` に追加するようにした。
  - これにより、非結合の 1-2 除外(`count_nonb_excl`)も、追加された結合項を通じて自動的に正しく設定される(除外リストは結合項から作られるため)。
  - CHARMM・GROMACS には同じ追加処理を入れていないため、`rigid_bond = YES` と固定原子の併用は AMBER 以外ではエラーにしている(`src/spdyn/sp_setup_spdyn.fpp` で `ctrl_data%ene_info%forcefield /= ForcefieldAMBER` を確認)。

### 4.5 その他の制限チェック

- `src/spdyn/spdyn.fpp`: 実行モードが `GenesisREMD` 以外で `fixatm_*` が指定されていればエラー。
- `src/spdyn/sp_setup_spdyn.fpp`: 並列 I/O(`pio_check_ranked_file` が真)、FEP(`alchemy` 引数あり)のときにエラー。

## 5. 制限事項(まとめ)

- gREST を含む REMD 専用。
- 並列 I/O の再開ファイル名(`(...)` を含む形式)は未対応。
- FEP(`[ALCHEMY]`)は未対応。
- 水分子は固定できない。
- `rigid_bond = YES` との併用は AMBER のみ。CHARMM・GROMACS は未対応。

## 6. 変更ファイル一覧

| ファイル | 内容 |
|---|---|
| `src/spdyn/sp_constraints_str.fpp` | `s_constraints` に `fixatm(:)` フラグ配列を追加 |
| `src/spdyn/sp_constraints.fpp` | `fixatm_select_index` / `fixatm_atom_index` の読み込み、`setup_fixatm_spdyn`、`apply_fixatm_domain`、`fixatm_requested` を追加 |
| `src/spdyn/sp_domain.fpp` | `unconstrain_fixed_hydrogen` を追加し、`setup_domain` の `check_light_atom_name` 直後から呼ぶ |
| `src/spdyn/sp_enefunc_amber.fpp` | `setup_enefunc_bond_constraint` で、SHAKE グループに一致しない H-bond を通常結合として追加 |
| `src/spdyn/sp_assign_velocity.fpp` | `stop_trans_rotation` で質量 0(固定)原子をスキップ |
| `src/spdyn/sp_setup_spdyn.fpp` | `setup_spdyn_remd` で固定原子の選択・適用・自由度減算・各種制限チェックを実施 |
| `src/spdyn/spdyn.fpp` | 実行モードによる制限チェック(REMD 以外は拒否) |

## 7. 検証

AMBER 版 `spdyn`(`--enable-gpu` なし)でビルドし、以下を確認した。

- CHARMM の DNA 系(T-REMD、交換なし/あり)で、固定原子(DNA の一部、94 原子)の変位が全フレームで 0、他の原子は動くこと。自由度が `3 × 固定原子数` 減ること。
- AMBER の溶媒和アラニンジペプチド系(`rigid_bond = YES`、固定原子 1-3)で、
  - 固定原子を含むグループの H-bond(例: 原子 2-4)が SHAKE から外れ、長さが変化すること。
  - 他の SHAKE された H-bond(例: 1.09 Å、0.9572 Å)は長さが一定のままであること。
  - 固定原子自体の変位が 0 であること。
  - SHAKE 本数が、固定原子を含むグループの分だけ正しく減ること。
- 実際の gREST 本番系(4YUS、151,113 原子、4 レプリカ、固定原子 106 個、剛体領域の固定)で、6000 ステップ実行後、固定領域の RMSD(アライン有無とも)が実質 0 であることを確認した(`ims` 上、`run_cpu` ディレクトリに検証スクリプト `check_fixatm_rmsd.py` と結果を保存済み)。
  - このとき `num_deg_freedom` のログで `(After setup of fixed atoms)` の行が `3 × 固定原子数` だけ減っていることも確認した。

開発初期に、固定原子を含む H-グループをそのまま `setup_hbond_group` で除外する方式(§4.4 とは別のアプローチ)を試したが、その方式ではドメイン振り分け(`sp_domain.fpp` の原子-セル割り当て)が H-グループ単位で動くため、グループから外れた原子がどのセルにも割り当てられず(`id_g2l` が 0 のまま)、計算から静かに消えるという不具合が見つかった。そのため、その方式は取り消し、§4.4 の「水素を軽い原子のフラグから外す」方式に切り替えている。

## 8. 今後の課題

- CHARMM・GROMACS での `rigid_bond = YES` と固定原子の併用(§4.4 と同様の結合項追加処理が必要)。
- 並列 I/O(PIO)経路への対応(再開ファイルからの質量・速度・自由度の設定)。
- FEP との併用。
- 長時間・交換ありでの安定性の確認(これまでの検証は数十 ps 〜 数十 ns 規模の短い実行)。
