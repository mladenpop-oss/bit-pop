# Bioconda Submit Instructions

## 1. Fork i clone bioconda-recipes

```bash
# Fork https://github.com/bioconda/bioconda-recipes na tvoj GitHub account
# Onda:

git clone https://github.com/mladenpop-oss/bioconda-recipes.git
cd bioconda-recipes
```

## 2. Kreiraj recipe folder

```bash
mkdir -p recipes/bit-pop
cp -r ~/path/to/bit-pop/bioconda-recipe/* recipes/bit-pop/
```

## 3. Provjeri strukturu

```
recipes/bit-pop/
├── meta.yaml
├── build.sh
└── bld.bat
```

## 4. Testiraj lokalno (opcionalno, preporučeno)

```bash
# Instaliraj conda-smithy (ako nemaš)
conda install -c conda-forge conda-smithy

# Testiraj recipe
conda build recipes/bit-pop

# Ili s proot za multi-platform:
conda install -c conda-forge anaconda-client proot
conda-build recipes/bit-pop --channel defaults --channel conda-forge --channel bioconda
```

## 5. Commit i push

```bash
git add recipes/bit-pop
git commit -m "add bit-pop recipe v0.2.0"
git push origin main
```

## 6. Kreiraj Pull Request

1. Idi na https://github.com/bioconda/bioconda-recipes
2. Klikni "Compare & pull request" (pojavit će se notifikacija za tvoj fork)
3. **PR title**: `bit-pop 0.2.0`
4. **PR description**:
```
Add recipe for bit-pop v0.2.0

Bit-Pop: Ultra-fast multi-genome DNA read classification
- FM-index + bit-level parallelism in Rust
- CAMI benchmark: 92.29% accuracy, 70% mapping
- Species-level: ~100% accuracy

Homepage: https://github.com/mladenpop-oss/bit-pop
Docker: ghcr.io/mladenpop-oss/bit-pop:latest
```

5. Klikni "Create pull request"

## 7. Čekaj CI review

Bioconda CI će automatski:
- Testirati build na Linux x86_64
- Provjeriti da li `bit-pop --help` radi
- Validirati meta.yaml format

Ako CI padne:
- Provjeri logs na GitHub Actions
- Fixaj i pushaj novi commit na istu branchu

## 8. Nakon merge-a

Kad maintainers mergeaju PR:
```bash
conda install -c bioconda bit-pop
```

Radi za sve!

## Važne napomene

- **Ne mijenjaj** `bioconda-recipe/` folder u originalnom bit-pop repo-u nakon submit-a
- Za novi release:
  1. Ažuriraj `version` i `sha256` u `meta.yaml`
  2. Ažuriraj `build: number` ako mijenjaš recipe bez nove verzije
  3. Kreiraj novi PR na bioconda-recipes
- Bioconda fokusira na **Linux x86_64** i **Linux aarch64**
- macOS/Windows support ide preko conda-forge ili Homebrew/Docker

## Linkovi

- Bioconda docs: https://bioconda.github.io/
- bioconda-recipes: https://github.com/bioconda/bioconda-recipes
- Recipe linter: https://bioconda.github.io/recipes.html
