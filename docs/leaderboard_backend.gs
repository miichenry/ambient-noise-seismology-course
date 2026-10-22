// ============================================================
// Google Apps Script — Leaderboard Backend
// ============================================================
// 1. Go to https://sheets.google.com → create a new blank sheet
// 2. Name it "ANT Quiz Leaderboard"
// 3. In the sheet, go to Extensions → Apps Script
// 4. Delete the default code and paste this entire file
// 5. Click Deploy → New deployment
//    - Type: Web app
//    - Execute as: Me
//    - Who has access: Anyone
// 6. Click Deploy → copy the URL
// 7. Paste the URL into quiz.html (replace BACKEND_URL)
// ============================================================

function doGet(e) {
  var sheet = SpreadsheetApp.getActiveSpreadsheet().getActiveSheet();
  var data = sheet.getDataRange().getValues();

  // Skip header row
  var results = [];
  for (var i = 1; i < data.length; i++) {
    if (data[i][0]) {
      results.push({
        name: data[i][0],
        score: data[i][1],
        total: data[i][2],
        pct: data[i][3],
        time: data[i][4],
        date: data[i][5]
      });
    }
  }

  // Sort by percentage desc, then time asc
  results.sort(function(a, b) {
    if (b.pct !== a.pct) return b.pct - a.pct;
    return a.time - b.time;
  });

  return ContentService
    .createTextOutput(JSON.stringify(results))
    .setMimeType(ContentService.MimeType.JSON);
}

function doPost(e) {
  var sheet = SpreadsheetApp.getActiveSpreadsheet().getActiveSheet();

  // Add header if sheet is empty
  if (sheet.getLastRow() === 0) {
    sheet.appendRow(['Name', 'Score', 'Total', 'Percentage', 'Time (s)', 'Date']);
  }

  var data = JSON.parse(e.postData.contents);

  sheet.appendRow([
    data.name,
    data.score,
    data.total,
    data.pct,
    data.time,
    new Date().toISOString()
  ]);

  return ContentService
    .createTextOutput(JSON.stringify({status: 'ok'}))
    .setMimeType(ContentService.MimeType.JSON);
}
