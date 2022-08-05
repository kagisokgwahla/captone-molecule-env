from .context import action
import unittest

class TestAction(unittest.TestCase):

    def test_setAction(self):
        act = action.Action()

        self.assertEqual(act.action_c, "")
        self.assertEqual(act.pos, "")
        self.assertEqual(act.mol, "")
        self.assertEqual(act.query, "")
        self.assertFalse(act.isSmarts, "")

        act.setAction("add", pos="front", mol="C")

        self.assertEqual(act.action_c, "add")
        self.assertEqual(act.pos, "front")
        self.assertEqual(act.mol, "C")
        self.assertEqual(act.query, "")
        self.assertFalse(act.isSmarts)
