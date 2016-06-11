define([
    'app',
    'text!app/templates/ui/tabs.hbs'
], function (App, template) {

    App.View.Ui = App.View.Ui || {};

    App.View.Ui.Tabs = App.Backbone.View.extend({
        tagName: 'div',
        className: 'list-tabs',

        template: App.Handlebars.compile(template),

        tabs: [],

        events: {
            "click a": "onTabClick"
        },

        initialize: function(options) {
            this.tabs = options.tabs || [];

            App._.each(this.tabs, function(tab, index) {
                tab.active = index == 0;
            });
        },

        setTabCount: function(name, count, add) {
            add = add || false;
            var tab;
            if (name) {
                tab = this._findTabByName(name);
            } else {
                tab = this.getActiveTab();
            }
            if (tab) {
                tab.count = add ? tab.count + count : count;
            }
            this.render();
        },

        getContext: function() {
            return this.getActiveTab().data;
        },

        getActiveTab: function() {
            return App._.find(this.tabs, function(tab) {
                return tab.active;
            });
        },

        _findTabByName: function(name) {
            return App._.find(this.tabs, function(tab) {
                return tab.name === name;
            });
        },

        onTabClick: function(e) {
            e.preventDefault();

            var tab = this._findTabByName(App.$(e.target).closest('a').data('name'));
            if (tab) {
                App._.each(this.tabs, function(tab) {
                    tab.active = false;
                });
                tab.active = true;
                this.render();
                this.trigger('context-change', tab.data);
            } else {
                throw 'Tab not found';
            }
        },

        render: function () {

            this.$el.html(this.template({tabs: this.tabs}));

            return this;
        }
    });

    App.View.Ui.Tabs.create = function(options, createOptions) {
        options = options || {};
        if (!options.tabs) {
            throw 'options.tabs is expected to be a hash of tabs';
        }

        var view = new App.View.Ui.Tabs(options);
        if (createOptions) {
            App.View.create(view, createOptions);
        }

        return view;
    };

    return App.View.Ui.Tabs;
});
