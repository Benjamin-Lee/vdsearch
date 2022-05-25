import logging

import rich_click as click
import typer
from click_didyoumean import DYMGroup
from rich_click import RichCommand as Command
from rich_click import RichGroup as Group

click.rich_click.USE_MARKDOWN = True
click.rich_click.MAX_WIDTH = 88

click.rich_click.COMMAND_GROUPS = {
    "vdsearch": [
        {"name": "Easy Workflows ([green]recommended[/])", "commands": ["easy-search"]},
        {
            "name": "circRNA detection and handling ([yellow]advanced[/])",
            "commands": [
                "find-circs",
                "canonicalize",
                "dedup",
            ],
        },
        {
            "name": "Ribozyme detection ([yellow]advanced[/])",
            "commands": ["infernal", "rnamotif", "ribozyme-filter"],
        },
        {
            "name": "Analysis Steps ([yellow]advanced[/])",
            "commands": [
                "fold",
                "orfs",
                "summarize",
                "cluster",
            ],
        },
        {
            "name": "Dataset Management",
            "commands": ["download-cms", "download-viroiddb", "purge"],
        },
        {"name": "Internal ([red]dangerous[/])", "commands": ["internal"]},
    ]
}
click.rich_click.STYLE_HELPTEXT = ""
click.rich_click.ERRORS_EPILOGUE = (
    "[dim]If you have any questions or need help, "
    "please contact me at [blue]https://benjamindlee.com\n"
)
click.rich_click.STYLE_FOOTER_TEXT = "dim"
click.rich_click.FOOTER_TEXT = (
    "For more information, please visit <https://benjamindlee.com>\n"
)


class MyRichGroup(Group, DYMGroup):
    """Wrapper to enable Did You Mean"""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


class MyTyper(typer.Typer):
    """A custom subclassed version of typer.Typer to allow rich help"""

    def __init__(self, *args, cls=MyRichGroup, **kwargs) -> None:
        super().__init__(*args, cls=cls, **kwargs)

    def command(self, *args, cls=Command, **kwargs) -> click.command:
        return super().command(*args, cls=cls, **kwargs)


def addLoggingLevel(levelName, levelNum, methodName=None):
    """
    Comprehensively adds a new logging level to the `logging` module and the
    currently configured logging class.

    `levelName` becomes an attribute of the `logging` module with the value
    `levelNum`. `methodName` becomes a convenience method for both `logging`
    itself and the class returned by `logging.getLoggerClass()` (usually just
    `logging.Logger`). If `methodName` is not specified, `levelName.lower()` is
    used.

    To avoid accidental clobberings of existing attributes, this method will
    raise an `AttributeError` if the level name is already an attribute of the
    `logging` module or if the method name is already present

    Example
    -------
    >>> addLoggingLevel('TRACE', logging.DEBUG - 5)
    >>> logging.getLogger(__name__).setLevel("TRACE")
    >>> logging.getLogger(__name__).trace('that worked')
    >>> logging.trace('so did this')
    >>> logging.TRACE
    5

    """
    if not methodName:
        methodName = levelName.lower()

    if hasattr(logging, levelName):
        raise AttributeError("{} already defined in logging module".format(levelName))
    if hasattr(logging, methodName):
        raise AttributeError("{} already defined in logging module".format(methodName))
    if hasattr(logging.getLoggerClass(), methodName):
        raise AttributeError("{} already defined in logger class".format(methodName))

    # This method was inspired by the answers to Stack Overflow post
    # http://stackoverflow.com/q/2183233/2988730, especially
    # http://stackoverflow.com/a/13638084/2988730
    def logForLevel(self, message, *args, **kwargs):
        if self.isEnabledFor(levelNum):
            self._log(levelNum, message, args, **kwargs)

    def logToRoot(message, *args, **kwargs):
        logging.log(levelNum, message, *args, **kwargs)

    logging.addLevelName(levelNum, levelName)
    setattr(logging, levelName, levelNum)
    setattr(logging.getLoggerClass(), methodName, logForLevel)
    setattr(logging, methodName, logToRoot)


addLoggingLevel("DONE", logging.INFO + 5)
