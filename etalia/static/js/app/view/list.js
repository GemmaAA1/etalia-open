define([
    'app',
    'text!app/templates/list.html'
], function (App, template) {

    //var defaults = {};


    return App.View.List = App.Backbone.View.extend({
        tagName: 'div',
        class: 'thumb-list',

        template: App.Handlebars.compile(template),

        events: {
            "click #thread-next-page": "onNextPageClick"
        },

        /*initialize: function (options) {
            //App.defaults(options, defaults);
        },*/

        render: function () {
            this.$el.html(this.template({}));

            this.mode.view
                .render()
                .appendTo(this.$('.document .inner'));

            return this;
        }
    });
});
