# Run the allocation test under a job object with a per-process memory cap, so that
# the allocation failure is real. Usage: powershell -File ... <cap in MB>
param([int]$capMb = 1200)
Add-Type -TypeDefinition @'
using System;
using System.Runtime.InteropServices;
public static class RevJob {
  [StructLayout(LayoutKind.Sequential)] public struct BASIC {
    public long PerProcessUserTimeLimit; public long PerJobUserTimeLimit;
    public uint LimitFlags; public UIntPtr MinimumWorkingSetSize; public UIntPtr MaximumWorkingSetSize;
    public uint ActiveProcessLimit; public UIntPtr Affinity; public uint PriorityClass; public uint SchedulingClass; }
  [StructLayout(LayoutKind.Sequential)] public struct IO {
    public ulong a, b, c, d, e, f; }
  [StructLayout(LayoutKind.Sequential)] public struct EXT {
    public BASIC Basic; public IO Io; public UIntPtr ProcessMemoryLimit; public UIntPtr JobMemoryLimit;
    public UIntPtr PeakProcessMemoryUsed; public UIntPtr PeakJobMemoryUsed; }
  [DllImport("kernel32.dll", CharSet = CharSet.Unicode)] public static extern IntPtr CreateJobObject(IntPtr a, string n);
  [DllImport("kernel32.dll")] public static extern bool SetInformationJobObject(IntPtr h, int c, ref EXT i, int l);
  [DllImport("kernel32.dll")] public static extern bool AssignProcessToJobObject(IntPtr h, IntPtr p);
  [DllImport("kernel32.dll")] public static extern IntPtr GetCurrentProcess();
  public static bool Cap(long bytes) {
    IntPtr h = CreateJobObject(IntPtr.Zero, null);
    EXT e = new EXT();
    e.Basic.LimitFlags = 0x100; // JOB_OBJECT_LIMIT_PROCESS_MEMORY
    e.ProcessMemoryLimit = new UIntPtr((ulong)bytes);
    if (!SetInformationJobObject(h, 9, ref e, Marshal.SizeOf(typeof(EXT)))) return false;
    return AssignProcessToJobObject(h, GetCurrentProcess());
  }
}
'@
$ok = [RevJob]::Cap([long]$capMb * 1MB)
"job cap $capMb MB assigned: $ok"
$env:PATH = "C:\rtools45\usr\bin;" + $env:PATH
& "C:\Program Files\R\R-4.6.1\bin\Rscript.exe" "C:\Users\adf44\source\r\frmtmb-wt-impsamp\dev\review-impsamp-fr1-alloc.R" 2>&1 | ForEach-Object { "$_" } | Where-Object { $_ -notmatch '^Boundary' }
"exit code: $LASTEXITCODE"
