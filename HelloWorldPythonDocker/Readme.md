1. 構建 Docker 鏡像

在終端（Terminal）中進入項目根目錄，執行以下命令來構建鏡像：

	docker build -t helloworld .

2. 運行容器

鏡像構建完成後，使用以下命令啟動容器：

	docker run --rm helloworld
	(--rm 參數表示容器運行結束後自動刪除，保持你的系統整潔。)

