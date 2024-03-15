del *.nupkg

nuget pack LinearAlgebra.csproj -IncludeReferencedProjects -Symbols


for %%f in (.\*.nupkg) do (
echo %%f|find ".symbols." >nul
if errorlevel 1 (echo %%f) else (nuget push -Source "https://ei-buildsrv.electroimpact.com/Automation/_packaging/AFP/nuget/v3/index.json" -ApiKey AzureDevOps %%f)
)

pause
