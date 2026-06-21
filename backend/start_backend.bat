@echo off
echo ========================================
echo  HEA Analysis Backend - Starting...
echo ========================================
cd /d "%~dp0"
uvicorn app.main:app --host 0.0.0.0 --port 8000 --reload
