# Reviewer of lane optrap, item 3: test-id-kron.R's merged fit and the
# test file, pinned to P cores (C03C03) and to E cores (3FC3FC).
# "none" as the library is a directory that does not exist, which
# .libPaths() drops, so that arm is rellib-r7.
Set-Location C:\Users\adf44\source\r\frmtmb-wt-optrap
$L = "C:/Users/adf44/source/r/wt-optrap-lib"
function Pin($mask, $name, $rest) {
  Start-Process -NoNewWindow -PassThru cmd.exe -ArgumentList "/c dev\optrap-rev-pin.cmd $mask $name $rest"
}
foreach ($i in 1..3) {
  $p = @()
  $p += Pin C03C03 "idk-laneP$i" "dev/optrap-idk.R $L laneP$i"
  $p += Pin 3FC3FC "idk-laneE$i" "dev/optrap-idk.R $L laneE$i"
  $p += Pin C03C03 "idk-r7P$i" "dev/optrap-idk.R none r7P$i"
  $p += Pin 3FC3FC "idk-r7E$i" "dev/optrap-idk.R none r7E$i"
  $p | Wait-Process
}
foreach ($i in 1..3) {
  $p = @()
  $p += Pin C03C03 "idkt-laneP$i" "dev/optrap-testfile.R $L frmtmb tests/testthat/test-id-kron.R"
  $p += Pin 3FC3FC "idkt-laneE$i" "dev/optrap-testfile.R $L frmtmb tests/testthat/test-id-kron.R"
  $p += Pin C03C03 "idkt-r7P$i" "dev/optrap-testfile.R none frmtmb tests/testthat/test-id-kron.R"
  $p += Pin 3FC3FC "idkt-r7E$i" "dev/optrap-testfile.R none frmtmb tests/testthat/test-id-kron.R"
  $p | Wait-Process
}
"IDKDONE"
