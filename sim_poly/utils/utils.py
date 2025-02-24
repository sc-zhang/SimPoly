import time


def time_print(info):
    print(
        "\033[32m%s\033[0m %s"
        % (time.strftime("[%H:%M:%S]", time.localtime(time.time())), info)
    )
