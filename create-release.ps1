<#
.SYNOPSIS
    Bumps the app version, commits, tags, and pushes - triggering the GitHub
    Actions release workflow (.github/workflows/release.yml) which builds the
    onefile exe and publishes the GitHub Release.

.EXAMPLE
    ./create-release.ps1 -Version 0.2.0
#>
param(
    [Parameter(Mandatory = $true)]
    [Alias("v")]
    [string]$Version
)

$ErrorActionPreference = "Stop"

function Fail($message) {
    Write-Error $message
    exit 1
}

if ($Version -notmatch '^\d+\.\d+\.\d+$') {
    Fail "Version must be in x.y.z form (got '$Version')."
}

if (-not (Get-Command git -ErrorAction SilentlyContinue)) {
    Fail "git not found on PATH."
}
if (-not (Get-Command uv -ErrorAction SilentlyContinue)) {
    Fail "uv not found on PATH."
}

$repoRoot = git rev-parse --show-toplevel 2>$null
if (-not $repoRoot) {
    Fail "Not inside a git repository."
}
Set-Location $repoRoot

$branch = git rev-parse --abbrev-ref HEAD
if ($branch -ne "main") {
    Fail "Must be on 'main' to release (currently on '$branch')."
}

$status = git status --porcelain
if ($status) {
    Fail "Working tree is not clean. Commit or stash changes first:`n$status"
}

$tag = "v$Version"
$existingTag = git tag --list $tag
if ($existingTag) {
    Fail "Tag '$tag' already exists."
}

Write-Host "Syncing with origin/main..."
git fetch origin
git pull --rebase origin main
if ($LASTEXITCODE -ne 0) { Fail "git pull --rebase failed. Resolve conflicts and re-run." }

Write-Host "Bumping version to $Version..."

$pyprojectPath = "pyproject.toml"
$pyproject = Get-Content $pyprojectPath -Raw
$updatedPyproject = $pyproject -replace '(?m)^version = "[^"]*"', "version = `"$Version`""
if ($updatedPyproject -eq $pyproject) {
    Fail "Could not find a 'version = \"...\"' line in $pyprojectPath."
}
Set-Content -Path $pyprojectPath -Value $updatedPyproject -NoNewline

$versionFilePath = "src/_version.py"
$versionFileContent = @"
"""Single source of truth for the app version, kept in sync with pyproject.toml
by create-release.ps1. Do not edit by hand outside that script.
"""

__version__ = "$Version"
"@
Set-Content -Path $versionFilePath -Value $versionFileContent -NoNewline

Write-Host "Syncing lockfile..."
uv lock
if ($LASTEXITCODE -ne 0) { Fail "uv lock failed." }

git add $pyprojectPath $versionFilePath uv.lock
git commit -m "chore: release $tag"
if ($LASTEXITCODE -ne 0) { Fail "git commit failed." }

Write-Host "Pushing main..."
git push origin main
if ($LASTEXITCODE -ne 0) { Fail "git push origin main failed." }

Write-Host "Tagging $tag..."
git tag $tag
git push origin $tag
if ($LASTEXITCODE -ne 0) { Fail "git push origin $tag failed." }

$remoteUrl = git config --get remote.origin.url
$slug = $remoteUrl -replace '^.*github\.com[:/]', '' -replace '\.git$', ''

Write-Host ""
Write-Host "Pushed $tag. GitHub Actions will build the exe and publish the release:"
Write-Host "  https://github.com/$slug/actions"
Write-Host "  https://github.com/$slug/releases/tag/$tag"
