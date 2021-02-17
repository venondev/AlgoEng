using Logging

function setupLogger()
    logger = SimpleLogger(stdout, Logging.Info)
    global_logger(logger)
end