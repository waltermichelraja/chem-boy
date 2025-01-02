import kivy
import sys
sys.path.append('../')
from kivy.app import App
from kivymd.theming import ThemeManager
from kivymd.uix.card import MDSeparator
from kivymd.uix.toolbar import MDToolbar
from kivy.uix.screenmanager import ScreenManager, Screen, NoTransition
from kivy.lang import Builder
from kivy.uix.floatlayout import FloatLayout
from kivy.uix.boxlayout import BoxLayout
from kivymd.uix.label import MDLabel
import Equations_V3 as eq
from kivy.properties import StringProperty
from kivy.core.text.markup import MarkupLabel
from kivy.uix.button import Button

# Fe2(SO4)3 + Ba(NO3)2 -> BaSO4 + Fe(NO3)3
# text_size: self.width, None


class MainMenu(Screen):
    pass


class WindowManager(ScreenManager):
    pass


class MainApp(App):
    theme_cls = ThemeManager()
    theme_cls.primary_palette = 'Pink'
    theme_cls.accent_palette = 'Pink'

    balanced_equation = StringProperty("")
    Button(text="Label is Added on screen !!:):)")

    def build(self):
        Button(text="Label is Added on screen !!:):)")
        return Builder.load_file("FinalUI.kv")

    def setEquation(self, etext):
        while etext.find("  ") > -1:
            etext = etext.replace("  ", " ")
        print(etext)

        fr = open("Equations.txt", "r+")
        if etext not in fr.read():
            fw = open("Equations.txt", "a+")
            fw.write(etext + "\r\n")
            fw.close()
        fr.close()

        return eq.balance(etext)


MainApp().run()