# 🔧 Complete Pipeline Troubleshooting Guide

## Issue: Complete Pipeline Feature Not Working

The backend API is confirmed working (tested successfully), so the issue is likely in the frontend.

## Quick Diagnosis Steps

### 1. Test the API Directly ✅
I've created a test page: `test_pipeline.html`

**To use it:**
1. Open `test_pipeline.html` in your browser
2. Click "Test Backend Health" - should show status "healthy"
3. Click "Test Complete Pipeline" - should show successful results
4. If this works, the issue is in the main frontend

### 2. Check Browser Console
1. Open your main app: http://localhost:3000
2. Press F12 → Console tab
3. Try running the complete pipeline
4. Look for specific error messages

### 3. Common Issues & Solutions

#### Issue A: Browser Cache
**Symptoms:** Old JavaScript, features not updating
**Solution:**
```
Hard refresh: Ctrl + Shift + R
Or clear browser cache completely
```

#### Issue B: CORS Issues
**Symptoms:** "CORS policy" errors in console
**Solution:** Backend already has CORS enabled, but try:
```
- Use http://localhost:3000 (not 127.0.0.1)
- Check if antivirus/firewall is blocking
```

#### Issue C: Network Connectivity
**Symptoms:** "Failed to fetch", "Network error"
**Solution:**
```
1. Verify backend is running: http://localhost:8000/docs
2. Verify frontend is running: http://localhost:3000
3. Check if ports are blocked
```

#### Issue D: JavaScript Runtime Errors
**Symptoms:** Console shows JavaScript errors
**Solution:**
```
1. Check console for specific errors
2. Look for undefined variables/functions
3. Check if React components are rendering
```

#### Issue E: Button Not Responding
**Symptoms:** Button clicks do nothing
**Solution:**
```
1. Check if protein is selected
2. Check if loading state is stuck
3. Verify onClick handlers are working
```

## Detailed Troubleshooting

### Step 1: Verify Backend
```powershell
# Test health
Invoke-WebRequest -Uri "http://localhost:8000/api/health"

# Test pipeline directly
$body = '{"sequence":"PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNF","name":"HIV-1 Protease","organism":"HIV-1","num_candidates":5}'
Invoke-WebRequest -Uri "http://localhost:8000/api/pipeline/complete" -Method POST -ContentType "application/json" -Body $body
```

### Step 2: Check Frontend State
In browser console, check:
```javascript
// Check if API constant is correct
console.log(UNIFIED_API); // Should be "http://localhost:8000/api"

// Check if protein is selected
console.log(selectedProtein); // Should not be empty

// Check current state
console.log(loading, progress, results);
```

### Step 3: Manual API Test
In browser console:
```javascript
// Test API call manually
fetch('http://localhost:8000/api/pipeline/complete', {
  method: 'POST',
  headers: { 'Content-Type': 'application/json' },
  body: JSON.stringify({
    sequence: "PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNF",
    name: "HIV-1 Protease",
    organism: "HIV-1",
    num_candidates: 5
  })
})
.then(r => r.json())
.then(data => console.log('✅ Success:', data))
.catch(err => console.error('❌ Error:', err));
```

## Expected Behavior

### When Working Correctly:
1. **Select protein** → Dropdown shows protein name
2. **Click "Run Complete Pipeline"** → Button shows loading spinner
3. **Progress bar** → Shows 0% → 10% → 50% → 100%
4. **Results appear** → Executive summary and detailed tabs
5. **Console shows:**
   ```
   ✅ Pipeline completed successfully
   Response data: {...}
   Has results? true
   ✅ State updated with results
   ```

### When Not Working:
- Button click does nothing
- Loading state gets stuck
- Error messages in console
- No results appear
- Alert with error message

## Quick Fixes to Try

### Fix 1: Restart Everything
```powershell
# Stop both servers (Ctrl+C in terminals)
# Then restart:
.\start-backend.ps1    # Terminal 1
.\start-frontend.ps1   # Terminal 2
```

### Fix 2: Clear Browser Data
```
1. Press Ctrl+Shift+Delete
2. Clear "Cached images and files"
3. Clear "Cookies and other site data"
4. Refresh page: Ctrl+Shift+R
```

### Fix 3: Try Different Browser
- Test in Chrome, Edge, or Firefox
- Use incognito/private mode
- Disable browser extensions

### Fix 4: Check Protein Selection
- Make sure a protein is selected from dropdown
- Try different example proteins
- Verify protein sequence is not empty

## Current System Status

✅ **Backend API:** Working (tested successfully)  
✅ **Endpoints:** All endpoints responding  
✅ **Data Format:** Correct response structure  
❓ **Frontend:** Needs diagnosis  

## Next Steps

1. **Use test_pipeline.html** to verify API works
2. **Check browser console** for specific errors
3. **Try hard refresh** (Ctrl+Shift+R)
4. **Test in different browser**
5. **Report specific error messages** you see

## Get Help

When reporting the issue, please include:
1. **Browser and version** (e.g., Chrome 119)
2. **Console error messages** (copy exact text)
3. **Steps you tried** from this guide
4. **What happens** when you click the button
5. **Results from test_pipeline.html**

---
**Created:** 2025-11-13 13:27 IST  
**Status:** Ready for diagnosis  
**Backend:** ✅ Confirmed working  
**Frontend:** 🔍 Needs investigation
