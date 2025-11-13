#!/bin/bash

# stAge macOS Build Script

set -e  # Exit on error

APP_NAME="stAge"
APP_PATH="dist/${APP_NAME}.app"
EXEC_PATH="${APP_PATH}/Contents/MacOS/${APP_NAME}"

echo "🔨 Building ${APP_NAME} for macOS..."
echo "🧹 Cleaning previous builds..."
rm -rf build dist *.spec

echo "📦 Running PyInstaller..."
pyinstaller --name="${APP_NAME}" \
  --windowed \
  --icon=icon.icns \
  --hidden-import=widgets \
  --hidden-import=themes \
  --hidden-import=plotting \
  --hidden-import=compute \
  --hidden-import=progress_dialog \
  --collect-submodules=scanpy \
  --collect-submodules=anndata \
  --copy-metadata=scikit-learn \
  --osx-entitlements-file=entitlements.plist \
  --clean \
  --noconfirm \
  main.py

# Check if build succeeded
if [ ! -d "${APP_PATH}" ]; then
  echo "❌ Build failed - ${APP_PATH} not found"
  exit 1
fi

echo "✅ Build completed"

echo "🔏 Code signing the app..."
codesign --force --deep --sign - \
  --entitlements entitlements.plist \
  --timestamp "${APP_PATH}"

echo "🔍 Verifying signature..."
codesign -dvvv "${APP_PATH}"

echo "🔍 Checking executable..."
if [ -f "${EXEC_PATH}" ]; then
  chmod +x "${EXEC_PATH}"
  ls -l "${EXEC_PATH}"
  echo "✅ Executable found and permissions set"
else
  echo "❌ Executable not found at ${EXEC_PATH}"
  ls -la "${APP_PATH}/Contents/MacOS/" || echo "MacOS directory doesn't exist"
  exit 1
fi

echo "🧪 Testing app launch (3s timeout)..."
EXIT_CODE=0

if command -v gtimeout &> /dev/null; then
  gtimeout 3s "${EXEC_PATH}" 2>&1 || EXIT_CODE=$?
else
  perl -e 'alarm shift; exec @ARGV' 3 "${EXEC_PATH}" 2>&1 || EXIT_CODE=$?
fi

if [ "$EXIT_CODE" -eq 142 ] || [ "$EXIT_CODE" -eq 124 ] || [ "$EXIT_CODE" -eq 14 ]; then
  echo "✅ App appears to launch successfully (timed out as expected)"
elif [ "$EXIT_CODE" -eq 0 ]; then
  echo "✅ App launched and exited normally"
else
  echo "⚠️  App exited with code $EXIT_CODE - check logs above"
fi

echo ""
echo "✅ Build complete: ${APP_PATH}"
echo ""
echo "To test manually:"
echo "  open ${APP_PATH}"
echo ""
echo "To view logs if app crashes:"
echo "  log show --predicate 'process == \"${APP_NAME}\"' --last 1m"
echo ""
echo "To check for missing dependencies:"
echo "  otool -L ${EXEC_PATH}"