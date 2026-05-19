# 使用轻量级 Python 镜像
FROM python:3.9-slim

# 设置容器内的工作目录
WORKDIR /app

# 将你的 app.py 复制到容器中
COPY hello.py .

# 运行程序
CMD ["python", "hello.py"]
