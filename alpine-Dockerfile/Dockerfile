# 使用官方的 GCC 映像檔作為基礎
FROM gcc:13

# 設定容器內的工作目錄
WORKDIR /usr/src/app

# 將當前資料夾的原始碼複製到容器中
COPY makefile *.c *.h ./

# 編譯 C 程式 (輸出成名為 main 的執行檔)
RUN make

# 容器啟動時執行編譯好的程式
CMD ["./main", "4", "0", "1"]
