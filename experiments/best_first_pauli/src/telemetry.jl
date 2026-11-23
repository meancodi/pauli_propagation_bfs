# External process monitoring; never equate resident bytes to allocation traffic.
function process_rss(pid::Integer)
    if Sys.iswindows()
        handle = ccall((:OpenProcess,"kernel32"),Ptr{Cvoid},(UInt32,Int32,UInt32),0x410,0,pid)
        handle == C_NULL && return missing
        data = zeros(UInt64,10)
        data[1] = 80
        ok = ccall((:GetProcessMemoryInfo,"psapi"),Int32,(Ptr{Cvoid},Ptr{UInt64},UInt32),handle,data,80)
        ccall((:CloseHandle,"kernel32"),Int32,(Ptr{Cvoid},),handle)
        return ok == 0 ? missing : Int(data[3])
    elseif Sys.islinux()
        try
            s = read("/proc/$pid/status",String)
            return parse(Int,match(r"VmRSS:\s+(\d+)",s).captures[1])*1024
        catch
            return missing
        end
    else
        try
            return parse(Int,strip(read(`ps -o rss= -p $pid`,String)))*1024
        catch
            return missing
        end
    end
end
