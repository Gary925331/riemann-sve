# ===== 建置階段 =====
FROM alpine:3.20 AS build
RUN apk add --no-cache build-base
WORKDIR /usr/src/app
COPY makefile *.c ./
RUN make

# ===== 執行階段 =====
FROM alpine:3.20
# OpenMP 執行期函式庫
RUN apk add --no-cache libgomp
WORKDIR /app
COPY --from=build /usr/src/app/main /app/main
# 用絕對路徑，working_dir 換掉也不會找不到
ENTRYPOINT ["/app/main"]
CMD ["1", "0", "1"]
