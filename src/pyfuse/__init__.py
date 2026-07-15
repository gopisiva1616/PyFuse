import logging
from datetime import datetime

__version__ = "0.1.0rc1"
class ConsoleFilter(logging.Filter):
    def filter(self, record: logging.LogRecord) -> bool:
        # default: show on console
        return getattr(record, "console", True)
    
def setup_logger(logger, log_path):

    log_name = '/pyfuse_' + datetime.now().strftime("%d-%m-%Y_%I-%M-%S") + '.log'
    log_file = log_path + log_name
    logger.setLevel(logging.DEBUG)

    console_handler = logging.StreamHandler()
    console_handler.setFormatter(logging.Formatter('%(asctime)s %(levelname)s [%(module)s]: %(message)s'))
    logger.addHandler(console_handler)
    console_handler.addFilter(ConsoleFilter())

    file_handler = logging.FileHandler(log_file)
    formatter = logging.Formatter('%(asctime)s %(levelname)s [%(module)s]: %(message)s')
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
